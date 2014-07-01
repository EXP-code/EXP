#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <map>
#include <algorithm>

#include "global.H"
#include "UserTreeDSMC.H"
#include "CollideIon.H"
#include "localmpi.h"
#include "Species.H"

using namespace std;

double   CollideIon::Nmin    = 1.0e-08;
double   CollideIon::Nmax    = 1.0e+25;
double   CollideIon::Tmin    = 1.0e+03;
double   CollideIon::Tmax    = 1.0e+08;
double   CollideIon::TolV    = 1.0e-03;
unsigned CollideIon::Nnum    = 400;
unsigned CollideIon::Tnum    = 200;
string   CollideIon::cache   = ".HeatCool";
bool     CollideIon::frost_warning = false; // For debugging . . . 

bool NO_COOL = false;

// Cross-section debugging, false for production
//
static bool CROSS_DBG  = false;

// Minimum energy for Rutherford scattering of ions used to estimate
// the elastic scattering cross section
//
const double FloorEv = 0.05;

CollideIon::CollideIon(ExternalForce *force, Component *comp, double hD, double sD, 
		       const std::string& smap, int Nth) : 
  Collide(force, comp, hD, sD, Nth)
{
  NUM = 0;

  parseSpecies(smap);

  // Fill the Chianti data base
  //
  ch.createIonList(ZList);

  // Cross-section storage
  //
  csections = std::vector<sKey2Dmap> (nthrds);

  // Random variable generators
  //
  gen  = new ACG(11+myid);
  unit = new Uniform(0.0, 1.0, gen);
  
  // Energy diagnostics
  //
  totalSoFar = 0.0;
  massSoFar  = 0.0;
  lostSoFar  = vector<double>(nthrds, 0.0);
  
  collD = boost::shared_ptr<collDiag>(new collDiag(this));

  if (myid==0 && NO_COOL) std::cout << "No cooling is set to TRUE" 
				    << std::endl;

  // Per thread workspace initialization
  //
  dCrossMap.resize(nthrds);
  dInterMap.resize(nthrds);
  sCrossMap.resize(nthrds);
  sInterMap.resize(nthrds);
  sWghtsMap.resize(nthrds);
  CE1      .resize(nthrds);
  CE2      .resize(nthrds);
  kCE1     .resize(nthrds);
  kCE2     .resize(nthrds);
  kEi      .resize(nthrds);
  kEe1     .resize(nthrds);
  kEe2     .resize(nthrds);
  Ein1     .resize(nthrds);
  Ein2     .resize(nthrds);
  spTau    .resize(nthrds);
  spCrm    .resize(nthrds);
  spProb   .resize(nthrds);

  //
  // Cross-section debugging [INIT]
  //
  if (CROSS_DBG) {
    nextTime_dbg = 0.0;		// Next target time
    nCnt_dbg     = 0;		// Number of cells accumulated so far

    ostringstream ostr;
    ostr << outdir << runtag << ".cross_section_dbg." << myid;
    cross_debug = ostr.str();
    std::ifstream in(cross_debug.c_str());
    if (!in) {
      std::ofstream out(cross_debug.c_str());
      out << std::setw( 8) << "Count"
	  << std::setw(18) << "Time"
	  << std::setw(18) << "Initial"
	  << std::setw(18) << "Final"
	  << std::setw(18) << "Ratio"
	  << std::endl
	  << std::setw( 8) << "-------"
	  << std::setw(18) << "-------"
	  << std::setw(18) << "-------"
	  << std::setw(18) << "-------"
	  << std::setw(18) << "-------"
	  << std::endl;
    }
  }

  // Enum collsion-type label fields
  //
  labels[geometric_1] = "geometric  [1]";
  labels[neut_elec_1] = "neutral el [1]";
  labels[ion_elec_1 ] = "charged el [1]";
  labels[free_free_1] = "free-free  [1]";
  labels[colexcite_1] = "col excite [1]";
  labels[ionize_1   ] = "ionization [1]";
  labels[recomb_1   ] = "recombine  [1]";
  labels[geometric_2] = "geometric  [2]";
  labels[neut_elec_2] = "neutral el [2]";
  labels[ion_elec_2 ] = "charged el [2]";
  labels[free_free_2] = "free-free  [2]";
  labels[colexcite_2] = "col excite [2]";
  labels[ionize_2   ] = "ionization [2]";
  labels[recomb_2   ] = "recombine  [2]";
}

CollideIon::~CollideIon()
{
}


/**
   Precompute all the necessary cross sections
 */
void CollideIon::initialize_cell(pHOT* tree, pCell* cell, double rvmax, int id)
{
				// Cache the calling tree
  curTree = tree;

  double KEtot, KEdspC;
  cell->KE(KEtot, KEdspC);	// KE in cell

  double massC = cell->Mass();	// Mass in cell

				// Used for diagnostics only
  totalSoFar += massC * KEdspC;
  massSoFar  += massC;
  
  // Representative avg cell velocity in cgs
  //  
  double vavg = 0.5*rvmax*UserTreeDSMC::Vunit;

  // Representative avg cell energy in ergs
  //
  double Eerg = 0.5*vavg*vavg*amu;

  // In eV
  //
  double EeV = Eerg / eV;

  // Upscaling factor for scattering cross section
  //
  double sUp  = diamfac*diamfac;

  if (aType == Direct) {

    typedef std::map<speciesKey, unsigned> Count;

    for (Count::iterator it1 = cell->count.begin(); it1 != cell->count.end(); it1++)  {

      speciesKey i1 = it1->first;
      double Cross1 = geometric(i1.first);
    
      for (Count::iterator 
	     it2 = cell->count.begin(); it2 != cell->count.end(); it2++) {
	
	speciesKey i2 = it2->first;
	double Cross2 = geometric(i2.first);

	double mu = atomic_weights[i1.first] * atomic_weights[i1.first] / 
	  (atomic_weights[i1.first] + atomic_weights[i2.first]);

	if (i2.second>1) {
	  double ne2   = i2.second - 1;
	  double eVel2 = sqrt(amu*atomic_weights[i2.second]/me);
	  if (i1.second==1)
	    Cross1 = elastic(i1.first, EeV * mu) * eVel2 * ne2;
	  else {
	    double b = 0.5*esu*esu*(i1.second - 1) /
	      std::max<double>(Eerg*mu, FloorEv*eV) * 1.0e7; // nm
	    Cross1 = M_PI*b*b * eVel2 * ne2;
	  }
	}

	if (i1.second>1) {
	  double ne1   = i1.second - 1;
	  double eVel1 = sqrt(amu*atomic_weights[i1.first]/me);
	  if (i2.second==1)
	    Cross2 = elastic(i2.first, EeV * mu) * eVel1 * ne1;
	  else {
	    double b = 0.5*esu*esu*(i2.second - 1) /
	      std::max<double>(Eerg*mu, FloorEv*eV) * 1.0e7; // nm
	    Cross2 = M_PI*b*b * eVel1 * ne1;
	  }
	}

	csections[id][i1][i2] = (Cross1 + Cross2) * sUp * 1e-14 / 
	  (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
      }
    }
  }

  if (aType == Trace) {

    csections[id][defaultKey][defaultKey] = 0.0;

    // Compute mean weights in the cell
    //
    typedef std::vector<unsigned long>::iterator bIter;
    typedef std::map<speciesKey, int>::iterator  sIter;

				// Mean weight
    keyWghtsMap meanW;
    for (sIter s=SpList.begin(); s!=SpList.end(); s++) 
      meanW[s->first] = sWghtsMap[id][s->first] = 0.0;

				// Total particle mass
    double massP = 0.0;
    double numbP = 0.0;
    
    for (bIter b=cell->bods.begin(); b!= cell->bods.end(); b++) {
				// Particle
      Particle *p = tree->Body(*b);
      massP += p->mass;
				// Mass-weighted trace fraction
      for (sIter s=SpList.begin(); s!=SpList.end(); s++) {
	meanW[s->first] += p->mass * p->dattrib[s->second];
	sWghtsMap[id][s->first] += p->dattrib[s->second]/atomic_weights[s->first.first];
      }
      numbP += 1.0;
    }
				// Compute mean weight
				//
    for (sIter s=SpList.begin(); s!=SpList.end(); s++) {
      meanW[s->first] /= massP;
      sWghtsMap[id][s->first] /= numbP;
    }

				// Compute cross sections
				//
    for (sIter s1=SpList.begin(); s1!=SpList.end(); s1++) {

      speciesKey i1 = s1->first;
      double Cross1 = geometric(i1.first);
      
      for (sIter s2=SpList.begin(); s2!=SpList.end(); s2++) {
	
	speciesKey i2 = s2->first;
	double Cross2 = geometric(i2.first);
	
	double mu = atomic_weights[i1.first] * atomic_weights[i1.first] / 
	  (atomic_weights[i1.first] + atomic_weights[i2.first]);

	if (i2.second>1) {
	  double ne2   = i2.second - 1;
	  double eVel2 = sqrt(amu*atomic_weights[i2.second]/me);
	  if (i1.second==1)
	    Cross1 = elastic(i1.first, EeV * mu) * eVel2 * ne2;
	  else {
	    double b = 0.5*esu*esu*(i1.second - 1) /
	      std::max<double>(Eerg*mu, FloorEv*eV) * 1.0e7; // nm
	    Cross1 = M_PI*b*b * eVel2 * ne2;
	  }
	}

	if (i1.second>1) {
	  double ne1   = i1.second - 1;
	  double eVel1 = sqrt(amu*atomic_weights[i1.first]/me);
	  if (i2.second==1)
	    Cross2 = elastic(i2.first, EeV * mu) * eVel1 * ne1;
	  else {
	    double b = 0.5*esu*esu*(i2.second - 1) /
	      std::max<double>(Eerg*mu, FloorEv*eV) * 1.0e7; // nm
	    Cross2 = M_PI*b*b * eVel1 * ne1;
	  }
	}

	double tCross = (Cross1 + Cross2) * 
	  meanW[i1] / atomic_weights[i1.first] *
	  meanW[i2] / atomic_weights[i2.first] *
	  sUp * 1e-14 / (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);

	csections[id][defaultKey][defaultKey] += tCross;
      }
    }
  }
}


sKey2Dmap& 
CollideIon::totalScatteringCrossSections(double crm, pCell *c, int id)
{
  typedef std::map<speciesKey, unsigned> Count;
  Count::iterator it1, it2;
  
  double vel  = crm * UserTreeDSMC::Vunit;
  double Eerg = 0.5*vel*vel*amu;
  double EeV  = Eerg / eV;

  // Upscaling factor for scattering cross section
  //
  double sUp  = diamfac*diamfac;
  
  if (aType == Direct) {

    for (it1 = c->count.begin(); it1 != c->count.end(); it1++)  {

      speciesKey i1 = it1->first;
      double geom1  = geometric(i1.first);
    
      for (it2 = c->count.begin(); it2 != c->count.end(); it2++)  
	{
	  speciesKey i2 = it2->first;
	  double geom2 = geometric(i2.first);
	  
	  double mu = atomic_weights[i1.first] * atomic_weights[i1.first] / 
	    (atomic_weights[i1.first] + atomic_weights[i2.first]);
	  
	  double eVel1 = sqrt(amu*atomic_weights[i1.first]/me);
	  double eVel2 = sqrt(amu*atomic_weights[i2.first]/me);
	  
	  double Cross1 = 0.0;
	  double Cross2 = 0.0;
	  
	  // Both particles neutral?
	  //
	  if (i1.second==1 and i2.second==2) {
	    Cross1 = geom1;
	    Cross2 = geom2;
	  }

	  // Electrons in second particle?
	  //
	  unsigned ne2 = i2.second - 1;
	  if (ne2) {
	    if (i1.second==1)	// Neutral atom-electron scattering
	      Cross1 = elastic(i1.first, EeV * mu) * eVel2*ne2;
	    else {		// Rutherford scattering
	      double b = 0.5*esu*esu*(i1.second - 1) /
		std::max<double>(Eerg*mu, FloorEv*eV) * 1.0e7; // nm
	      Cross1 = M_PI*b*b * eVel2*ne2;
	    }
	  }

	  // Electrons in first particle?
	  //
	  unsigned ne1 = i1.second - 1;
	  if (ne1) {
	    if (i2.second==1)	// Neutral atom-electron scattering
	      Cross2 = elastic(i2.first, EeV * mu) * eVel1*ne1;
	    else {		// Rutherford scattering
	      double b = 0.5*esu*esu*(i2.second - 1) /
		std::max<double>(Eerg*mu, FloorEv*eV) * 1.0e7; // nm
	      Cross2 = M_PI*b*b * eVel1*ne1;
	    }
	  }
	  
	  csections[id][i1][i2] = (Cross1 + Cross2) * sUp * 1e-14 / 
	    (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
	}
    }

  }

  
  if (aType == Trace) {

    csections[id][defaultKey][defaultKey] = 0.0;

    // Compute the mean trace weight in the cell
    //
    typedef std::vector<unsigned long>::iterator bIter;
    typedef std::map<speciesKey, int>::iterator  sIter;

    std::map<speciesKey, double> meanW;
    for (sIter s=SpList.begin(); s!=SpList.end(); s++)
      meanW[s->first] = sWghtsMap[id][s->first] = 0.0;

    double massP = 0.0;
    double numbP = 0.0;

    for (bIter b=c->bods.begin(); b!= c->bods.end(); b++) {
				// Particle
      Particle *p = curTree->Body(*b);
				// Mass accumulation
      massP += p->mass;
				// Mean weight accumulation
      for (sIter s=SpList.begin(); s!=SpList.end(); s++) {
	meanW[s->first] += p->mass * p->dattrib[s->second];
	sWghtsMap[id][s->first] += p->dattrib[s->second]/atomic_weights[s->first.first];
      }
      numbP += 1.0;
    }
    
				// Normalize
    for (sIter s=SpList.begin(); s!=SpList.end(); s++) {
      meanW[s->first] /= massP;
      sWghtsMap[id][s->first] /= numbP;
    }

				// Compute cross sections
				//
    for (sIter s1=SpList.begin(); s1!=SpList.end(); s1++) {

      speciesKey i1 = s1->first;
      double geom1  = geometric(i1.first);
    
      for (sIter s2=SpList.begin(); s2!=SpList.end(); s2++) {

	speciesKey i2 = s2->first;
	double geom2 = geometric(i2.first);
	  
	double mu = atomic_weights[i1.first] * atomic_weights[i1.first] / 
	  (atomic_weights[i1.first] + atomic_weights[i2.first]);
	  
	double eVel1 = sqrt(amu*atomic_weights[i1.first]/me);
	double eVel2 = sqrt(amu*atomic_weights[i2.first]/me);
	  
	double Cross1 = 0.0;
	double Cross2 = 0.0;
	  
	// Both particles neutral?
	//
	if (i1.second==1 and i2.second==1) {
	  Cross1 = geom1;
	  Cross2 = geom2;
	}
	
	// Electrons in second particle?
	//
	unsigned ne2 = i2.second - 1;
	if (ne2) {
	  if (i1.second==1)	// Neutral atom-electron scattering
	    Cross1 = elastic(i1.first, EeV * mu) * eVel2 * ne2;
	  else {		// Rutherford scattering
	    double b = 0.5*esu*esu*(i1.second - 1) /
	      std::max<double>(Eerg*mu, FloorEv*eV) * 1.0e7; // nm
	    Cross1 = M_PI*b*b * eVel2 * ne2;
	  }
	}

	// Electrons in first particle?
	//
	unsigned ne1 = i1.second - 1;
	if (ne1) {
	  if (i2.second==1)	// Neutral atom-electron scattering
	    Cross2 = elastic(i2.first, EeV * mu) * eVel1 * ne1;
	  else {		// Rutherford scattering
	    double b = 0.5*esu*esu*(i2.second - 1) /
	      std::max<double>(Eerg*mu, FloorEv*eV) * 1.0e7; // nm
	    Cross2 = M_PI*b*b * eVel1 * ne1;
	  }
	}
	
	double tCross = (Cross1 + Cross2) *
	  meanW[i1] / atomic_weights[i1.first] *
	  meanW[i2] / atomic_weights[i2.first] *
	  sUp * 1e-14 / (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);

	csections[id][defaultKey][defaultKey] += tCross;
      }
    }
  }

  return csections[id];
}

double CollideIon::crossSectionDirect(pHOT *tree, Particle* p1, Particle* p2, 
				      double cr, int id)
{
  // Species keys
  //
  KeyConvert k1(p1->iattrib[use_key]), k2(p2->iattrib[use_key]);

  unsigned short Z1 = k1.getKey().first, C1 = k1.getKey().second;
  unsigned short Z2 = k2.getKey().first, C2 = k2.getKey().second;
  
  // Number of atoms in each super particle
  //
  double N1 = (p1->mass*UserTreeDSMC::Munit)/(atomic_weights[Z1]*amu);
  double N2 = (p2->mass*UserTreeDSMC::Munit)/(atomic_weights[Z2]*amu);

  // Number of associated electrons for each particle
  //
  double ne1 = C1 - 1;
  double ne2 = C2 - 1;
  
  // Energy available in the center of mass of the atomic collision
  //
  double m1  = atomic_weights[Z1]*amu;
  double m2  = atomic_weights[Z2]*amu;
  double mu  = m1 * m2 / (m1 + m2);
  double vel = cr * UserTreeDSMC::Vunit;

  // Translational COM energy
  //
  kEi[id] = 0.5 * mu * vel*vel;

  // Electron velocity equipartition factors
  //
  double eVel1 = sqrt(m1/me);
  double eVel2 = sqrt(m2/me);

  // Internal energy per particle
  //
  Ein1[id] = Ein2[id] = 0.0;
  if (use_Eint>=0) {
    Ein1[id] = p1->dattrib[use_Eint] * UserTreeDSMC::Eunit/N1;
    Ein2[id] = p2->dattrib[use_Eint] * UserTreeDSMC::Eunit/N2;

    // Compute the total available energy and divide among degrees of freedom
    // Convert ergs to eV
    //
    kEe1[id] = (kEi[id] + Ein1[id])/(1.0 + ne1) / eV;
    kEe2[id] = (kEi[id] + Ein2[id])/(1.0 + ne2) / eV;
  } else {
    kEe1[id] = kEi[id] / eV;
    kEe2[id] = kEi[id] / eV;
  }
  
  
  /***
      INSERT HERE THE CALCULATIONS FOR DELTA E USING THE ION CROSS SECTIONS
  ***/

  /***
      Interaction integers:
      1: p1 ff
      2: p1 CE
      3: p1 DI
      4: p1 RE
      5: p2 ff
      6: p2 CE
      7: p2 DI
      8: p2 RE
  ***/
  
  // Save the per-interaction cross sections
  dCrossMap[id] = std::vector<double>();

  // Index the interactions
  dInterMap[id] = std::vector<int   >();

  double sum12 = 0.0;		// Accumulate inelastic total cross
  double sum21 = 0.0;		// sections as we go

  
  //--------------------------------------------------
  // Total scattering cross section
  //--------------------------------------------------

  double cross12 = 0.0;
  double cross21 = 0.0;

				//-------------------------------
				// Both particles neutral
				//-------------------------------
  if (C1==1 and C2==2) {
				// Geometric cross sections based on
				// atomic radius
    cross12 = geometric(Z1);
    dCrossMap[id].push_back(cross12*diamfac*diamfac);
    dInterMap[id].push_back(geometric_1);

    cross21 = geometric(Z2);
    dCrossMap[id].push_back(cross21*diamfac*diamfac);
    dInterMap[id].push_back(geometric_2);
  }

				//-------------------------------
				// Electrons in second particle
				//-------------------------------
  if (ne2 > 0) {
    if (C1==1) {		// Neutral atom-electron scattering
      cross12 = elastic(Z1, kEe2[id]) * eVel2*ne2;
      dCrossMap[id].push_back(cross12);
      dInterMap[id].push_back(neut_elec_1);
    }  else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C1-1) /
	std::max<double>(kEe2[id]*eV, FloorEv*eV) * 1.0e7; // nm
      cross12 = M_PI*b*b * eVel2*ne2;
      dCrossMap[id].push_back(cross12);
      dInterMap[id].push_back(ion_elec_1);
    }
  }
    
				//-------------------------------
				// Electrons in first particle
				//-------------------------------
  if (ne1 > 0) {
    if (C2==1) {		// Neutral atom-electron scattering
      cross21 = elastic(Z2, kEe1[id]) * eVel1*ne1;
      dCrossMap[id].push_back(cross21);
      dInterMap[id].push_back(neut_elec_2);
    } else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C2-1) /
	std::max<double>(kEe1[id]*eV, FloorEv*eV) * 1.0e7; // nm
      cross21 = M_PI*b*b * eVel1*ne1;
      dCrossMap[id].push_back(cross21);
      dInterMap[id].push_back(ion_elec_2);
    }
  }

  //--------------------------------------------------
  // Ion keys
  //--------------------------------------------------

  lQ Q1(Z1, C1), Q2(Z2, C2);

  //--------------------------------------------------
  // Particle 1 interacts with Particle 2
  //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
  if (C1 > 1 and ne2 > 0) {	// Ion and Ion only

    double ff1 = ch.IonList[Q1].freeFreeCross(kEe2[id]);
    double crs = eVel2*ne2 * ff1;

    if (crs>0.0) {
      dCrossMap[id].push_back(crs);
      dInterMap[id].push_back(free_free_1);
      sum12 += crs;
    }
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  if (ne2 > 0 and C1 <= Z1) {	// Particle 1 must be bound

    CE1[id] = ch.IonList[Q1].collExciteCross(kEe2[id]);

    double crs = eVel2*ne2 * CE1[id].back().first;

    if (crs>0.0) {
      dCrossMap[id].push_back(crs);
      dInterMap[id].push_back(colexcite_1);
      sum12 += crs;
    }
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  if (ne2 > 0 and C1 <= Z1) {	// Particle 1 must be bound

    double DI1 = ch.IonList[Q1].directIonCross(kEe2[id]);
    double crs = eVel2*ne2 * DI1;

    if (crs>0.0) {
      dCrossMap[id].push_back(crs);
      dInterMap[id].push_back(ionize_1);
      sum12 += crs;
    }
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  if (C1 > 1 and ne2 > 0) {	// Particle 1 must be an ion

    std::vector<double> RE1 = ch.IonList[Q1].radRecombCross(kEe2[id]);
    double crs = eVel2*ne2 * RE1.back();

    if (crs>0.0) {
      dCrossMap[id].push_back(crs);
      dInterMap[id].push_back(recomb_1);
      sum12 += crs;
    }
  }

  
  //--------------------------------------------------
  // Particle 2 interacts with Particle 1
  //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
  if (C2 > 1 and ne1 > 0) {
    double ff2 = ch.IonList[Q2].freeFreeCross(kEe1[id]);
    double crs = eVel1*ne1 * ff2;

    if (crs>0.0) {
      dCrossMap[id].push_back(crs);
      dInterMap[id].push_back(free_free_2);
      sum21 += crs;
    }
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  if (ne1 > 0 and C2 <= Z2) {

    CE2[id] = ch.IonList[Q2].collExciteCross(kEe1[id]);
    double crs = eVel1*ne1 * CE2[id].back().first;

    if (crs>0.0) {
      dCrossMap[id].push_back(crs);
      dInterMap[id].push_back(colexcite_2);
      sum21 += crs;
    }
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  if (ne1 > 0 and C2 <= Z2) {
    double DI2 = ch.IonList[Q2].directIonCross(kEe1[id]);
    double crs = ne1 * DI2;

    if (crs>0.0) {
      dCrossMap[id].push_back(crs);
      dInterMap[id].push_back(ionize_2);
      sum21 += crs;
    }
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  if (C2 > 1 and ne1 > 0) {
    std::vector<double> RE2 = ch.IonList[Q2].radRecombCross(kEe1[id]);
    double crs = eVel1*ne1*RE2.back();

    if (crs>0.0) {
      dCrossMap[id].push_back(crs);
      dInterMap[id].push_back(recomb_2);
      sum21 += crs;
    } 
  }

				//-------------------------------
				// *** Convert to system units
				//-------------------------------
  return (cross12 + cross21 + sum12 + sum21) * 1e-14 / 
    (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
}


double CollideIon::crossSectionTrace(pHOT *tree, Particle* p1, Particle* p2, 
				     double cr, int id)
{
  typedef std::map<speciesKey, int>::iterator sIter;

  double totalCross = 0.0;

  // Compute "matrix" of interactions
  //
  sIter sEnd(SpList.end());
  for (sIter s1=SpList.begin(); s1!=sEnd; s1++) {

				// Particle 1: key and number weight
    speciesKey k1 = s1->first;
    double     w1 = p1->dattrib[s1->second]/atomic_weights[k1.first];
    
    for (sIter s2=SpList.begin(); s2!=sEnd; s2++) {

				// Particle 1: key and number weight
      speciesKey k2 = s2->first;
      double     w2 = p2->dattrib[s2->second]/atomic_weights[k2.first];

				// Key of interaction pair
      dKey dkey(k1, k2);

				// Weight product
      double     ww = w1 * w2;

				// Atomic numbers
      unsigned short Z1 = k1.first, C1 = k1.second;
      unsigned short Z2 = k2.first, C2 = k2.second;
  

      // Number of associated electrons for each particle
      //
      double ne1 = C1 - 1;
      double ne2 = C2 - 1;
  
      // Energy available in the center of mass of the atomic collision
      //
      double m1  = amu;
      double m2  = amu;

      double mu  = m1 * m2 / (m1 + m2);
      double vel = cr * UserTreeDSMC::Vunit;

      // Translational COM energy
      //
      kEi[id] = 0.5 * mu * vel*vel;

      // Electron velocity equipartition factors
      //
      double eVel1 = sqrt(m1/me);
      double eVel2 = sqrt(m2/me);

      // Convert the total available energy from ergs to eV
      //
      
      kEe1[id] = kEi[id] / eV;
      kEe2[id] = kEi[id] / eV;
  
      /***
	  INSERT HERE THE CALCULATIONS FOR DELTA E USING THE ION CROSS SECTIONS
      ***/
      
      /***
	  Interaction integers:
	  1: p1 ff
	  2: p1 CE
	  3: p1 DI
	  4: p1 RE
	  5: p2 ff
	  6: p2 CE
	  7: p2 DI
	  8: p2 RE
      ***/
  
      // Save the per-interaction cross sections
      std::vector<double> tCrossMap;
      
      // Index the interactions
      std::vector<int> tInterMap;

      double sum12 = 0.0;	// Accumulate inelastic total cross
      double sum21 = 0.0;	// sections as we go

  
      //--------------------------------------------------
      // Total scattering cross section
      //--------------------------------------------------

      double cross12 = 0.0;
      double cross21 = 0.0;

      //-------------------------------
      // Elastic: both particles
      // are neutral
      //-------------------------------
      if (C1==1 and C2==1) {
	double sUp = diamfac * diamfac;
				// Geometric cross sections based on
				// atomic radius
	cross12 = geometric(Z1) * ww;
	tCrossMap.push_back(cross12 * sUp);
	tInterMap.push_back(geometric_1);

	cross21 = geometric(Z2) * ww;
	tCrossMap.push_back(cross21 * sUp);
	tInterMap.push_back(geometric_2);
      }

      //-------------------------------
      // Elastic: first particle is
      // neutral, electrons in second 
      // particle
      //-------------------------------
      if (ne2 > 0) {
				//-------------------------------
	if (C1==1) {		// *** Neutral atom-electron 
				// *** scattering
				//-------------------------------
	  cross12 = elastic(Z1, kEe2[id]) * eVel2 * ww * ne2;
	  tCrossMap.push_back(cross12);
	  tInterMap.push_back(neut_elec_1);

	}  else {			
				//-------------------------------
				// *** Rutherford scattering
				//-------------------------------
	  double b = 0.5*esu*esu*(C1-1) /
	    std::max<double>(kEe2[id]*eV, FloorEv*eV) * 1.0e7; // nm
	  cross12 = M_PI*b*b * eVel2 * ww * ne2;
	  tCrossMap.push_back(cross12);
	  tInterMap.push_back(ion_elec_1);
	}
      }
    
      //-------------------------------
      // Elastic: second particle is
      // neutral, electrons in first
      // particle
      //-------------------------------
      if (ne1 > 0) {
				//-------------------------------
	if (C2==1) {		// *** Neutral atom-electron 
				// *** scattering
				//-------------------------------
	  cross21 = elastic(Z2, kEe1[id]) * eVel1 * ww * ne1;
	  tCrossMap.push_back(cross21);
	  tInterMap.push_back(neut_elec_2);
	} else {
				//-------------------------------
				// *** Rutherford scattering
				//-------------------------------
	  double b = 0.5*esu*esu*(C2-1) /
	    std::max<double>(kEe1[id]*eV, FloorEv*eV) * 1.0e7; // nm
	  cross21 = M_PI*b*b * eVel1 * ww * ne1;
	  tCrossMap.push_back(cross21);
	  tInterMap.push_back(ion_elec_2);
	}
      }

      //--------------------------------------------------
      // Ion keys
      //--------------------------------------------------
      lQ Q1(Z1, C1), Q2(Z2, C2);

      //--------------------------------------------------
      // Inelastic scattering: Particle 1 interacts with 
      // electrons in Particle 2
      //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
      if (C1 > 1 and ne2 > 0) {	// Ion and Ion only
	double ff1 = ch.IonList[Q1].freeFreeCross(kEe2[id]);
	double crs = eVel2*ne2 * ff1 * ww;

	if (crs>0.0) {
	  tCrossMap.push_back(crs);
	  tInterMap.push_back(free_free_1);
	  sum12 += crs;
	}
      }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
      if (ne2 > 0 and C1 <= Z1) { // Particle 1 must be bound

	CEvector V = ch.IonList[Q1].collExciteCross(kEe2[id]);
	double crs = eVel2*ne2 * V.back().first * ww;

	if (crs>0.0) {
	  tCrossMap.push_back(crs);
	  tInterMap.push_back(colexcite_1);
	  kCE1[id][dkey] = V;
	  sum12 += crs;
	}
      }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
      if (ne2 > 0 and C1 <= Z1) { // Particle 1 must be bound

	double DI1 = ch.IonList[Q1].directIonCross(kEe2[id]);
	double crs = eVel2*ne2 * DI1 * ww;

	if (crs>0.0) {
	  tCrossMap.push_back(crs);
	  tInterMap.push_back(ionize_1);
	  sum12 += crs;
	}
      }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
      if (C1 > 1 and ne2 > 0) {	// Particle 1 must be an ion

	std::vector<double> RE1 = ch.IonList[Q1].radRecombCross(kEe2[id]);
	double crs = eVel2*ne2 * RE1.back() * ww;

	if (crs>0.0) {
	  tCrossMap.push_back(crs);
	  tInterMap.push_back(recomb_1);
	  sum12 += crs;
	}
      }

  
      //--------------------------------------------------
      // Inelastic scattering: Particle 2 interacts with 
      // electrons in Particle 1
      //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
      if (C2 > 1 and ne1 > 0) {
	double ff2 = ch.IonList[Q2].freeFreeCross(kEe1[id]);
	double crs = eVel1*ne1 * ff2 * ww;

	if (crs>0.0) {
	  tCrossMap.push_back(crs);
	  tInterMap.push_back(free_free_2);
	  sum21 += crs;
	}
      }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
      if (ne1 > 0 and C2 <= Z2) {

	CEvector V = ch.IonList[Q2].collExciteCross(kEe1[id]);
	double crs = eVel1*ne1 * V.back().first * ww;

	if (crs>0.0) {
	  tCrossMap.push_back(crs);
	  tInterMap.push_back(colexcite_2);
	  kCE2[id][dkey] = V;
	  sum21 += crs;
	}

      }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
      if (ne1 > 0 and C2 <= Z2) {
	double DI2 = ch.IonList[Q2].directIonCross(kEe1[id]);
	double crs = ne1 * DI2 * ww;
	
	if (crs>0.0) {
	  tCrossMap.push_back(crs);
	  tInterMap.push_back(ionize_2);
	  sum21 += crs;
	}
      }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
      if (C2 > 1 and ne1 > 0) {
	std::vector<double> RE2 = ch.IonList[Q2].radRecombCross(kEe1[id]);
	double crs = eVel1*ne1*RE2.back() * ww;

	if (crs>0.0) {
	  tCrossMap.push_back(crs);
	  tInterMap.push_back(recomb_2);
	  sum21 += crs;
	}
      } 

      sCrossMap[id][dkey] = tCrossMap;
      sInterMap[id][dkey] = tInterMap;

				//-------------------------------
				// *** Convert to system units
				//-------------------------------

      double tCross = (cross12 + cross21 + sum12 + sum21) * 1e-14 / 
	(UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
      
      totalCross += tCross;
    }
  }
      
  return totalCross;
}


int CollideIon::inelasticDirect(pHOT *tree, Particle* p1, Particle* p2, 
				double *cr, int id)
{
  int ret = 0;			// No error (flag)
  int interFlag = -1;

  // Species keys
  //
  KeyConvert k1(p1->iattrib[use_key]), k2(p2->iattrib[use_key]);
  collTDPtr ctd1 = (*collD)[k1.getKey()], ctd2 = (*collD)[k2.getKey()];
  unsigned short Z1 = k1.getKey().first, C1 = k1.getKey().second;
  unsigned short Z2 = k2.getKey().first, C2 = k2.getKey().second;
  
  // Number of atoms in each super particle
  //
  double N1 = (p1->mass*UserTreeDSMC::Munit)/(atomic_weights[Z1]*amu);
  double N2 = (p2->mass*UserTreeDSMC::Munit)/(atomic_weights[Z2]*amu);
  double NN = std::min<double>(N1, N2);	// Currently, N1 should equal N2
  
  // Number of associated electrons for each particle
  //
  double ne1 = C1 - 1;
  double ne2 = C2 - 1;

  // The total mass in system units
  //
  double Mt = p1->mass + p2->mass;
  if (Mt<=0.0) return ret;
  
  // Reduced mass in ballistic collision (system units)
  //
  double Mu = p1->mass * p2->mass / Mt;

  // Center of mass energy in the ballistic collision (system units)
  //
  double kE  = 0.5*Mu*(*cr)*(*cr);

  // For tracking energy conservation (system units)
  //
  double dE   = kE*TolV*TolV;
  double delE  = 0.0, delEeV = 0.0;

  // Now that the interactions have been calculated, create the
  // normalized cross section list to pick the interaction
  //
  std::vector<double> TotalCross;
  double tCross = 0.0;
  int si = 0;
  for (size_t i = 0; i < dCrossMap[id].size(); i++) {
    tCross += dCrossMap[id][i];
    TotalCross.push_back(tCross);
    si++;
  }

  assert (TotalCross.size() == dCrossMap[id].size());

  int partflag = 0;		// Will be 1 or 2, dependending on
				// which ion or neutral is selected
				// for interaction.  Will be 0 if no
				// interaction is selected.
  
  if (tCross != 0) {
    // Cumulative cross-section distribution for interaction selection
    //
    std::vector<double> CDF;
    for (size_t i = 0; i < TotalCross.size(); i++) {
      CDF.push_back(TotalCross[i]/tCross);
    }
    
    // Use a random variate to select the interaction from the
    // discrete cumulatative probability distribution (CDF)
    //
    double ran = (*unit)();
    int index  = -1;
    for (size_t i = 0; i < CDF.size(); i++) {
      if (ran < CDF[i]) {
	index = static_cast<int>(i);
	break;
      }
    }

    // Sanity check: did not assign index??
    //
    if (index<0) {
      std::cout << "CDF location falure, myid=" << myid
		<< ", ran=" << ran << std::endl;
      index = 0;
    }

    // Finally, set the interaction type based on the selected index
    //
    interFlag = dInterMap[id][index];

    //-------------------------
    // VERBOSE DEBUG TEST
    //-------------------------
    // Set to false for production
    //          |
    //          v
    static bool DEBUG_F = false;
    //
    if (DEBUG_F) {
      //
      // Output on collisions for now . . . 
      //
      if (interFlag % 100 == 4) {
	std::cout << std::setw( 8) << "index"
		  << std::setw( 8) << "flag"
		  << std::setw(14) << "cross"
		  << std::setw(14) << "cumul"
		  << std::setw(18) << "type label"
		  << std::endl
		  << std::setw( 8) << "-----"
		  << std::setw( 8) << "-----"
		  << std::setw(14) << "---------"
		  << std::setw(14) << "---------"
		  << std::setw(18) << "---------------"
		  << std::endl;
	for (size_t i = 0; i < dCrossMap[id].size(); i++) {
	  std::cout << std::setw( 8) << i
		    << std::setw( 8) << dInterMap[id][i]
		    << std::setw(14) << dCrossMap[id][i]
		    << std::setw(14) << CDF[i]
		    << std::setw(18) << labels[dInterMap[id][i]]
		    << std::endl;
	}
	std::cout << std::endl;
      }
    }

    //--------------------------------------------------
    // Ion keys
    //--------------------------------------------------

    lQ Q1(Z1, C1), Q2(Z2, C2);

    //-------------------------
    // Particle 1 interactions
    //-------------------------

    if (interFlag == free_free_1) {
      delE          = IS.selectFFInteract(ch.IonList[Q1], kEe2[id]);
      partflag      = 1;
      ctd1->ff[id].first++; 
      ctd1->ff[id].second  += delE;
    }

    if (interFlag == colexcite_1) {
      delE = IS.selectCEInteract(ch.IonList[Q1], CE1[id]);
      partflag      = 1;
      ctd1->CE[id].first++; 
      ctd1->CE[id].second  += delE;
    }

    if (interFlag == ionize_1) {
      delE          = IS.DIInterLoss(ch.IonList[Q1]);
      p1->iattrib[use_key] = k1.updateC(++C1);
      assert(C1 <= (Z1 + 1));
      partflag      = 1;
      ctd1->CI[id].first++; 
      ctd1->CI[id].second  += delE;
    }

    // KE carried by electron is subtracted from the thermal reservoir
    // but not the radiation of "binding".  The "binding" radiation
    // decreases the total energy of the gas but not the thermal
    // component.
    //
    if (interFlag == recomb_1) {
      delE          = kEe2[id];
      p1->iattrib[use_key] = k1.updateC(--C1);
      assert(C1 > 0);
      partflag      = 1;
      ctd1->RR[id].first++; 
      ctd1->RR[id].second  += delE;
    }
    
    //-------------------------
    // Particle 2 interactions
    //-------------------------

    if (interFlag == free_free_2) {
      delE          = IS.selectFFInteract(ch.IonList[Q2], kEe1[id]);
      partflag      = 2;
      ctd2->ff[id].first++;
      ctd2->ff[id].second += delE;
    }

    if (interFlag == colexcite_2) {
      delE         = IS.selectCEInteract(ch.IonList[Q2], CE2[id]);
      partflag     = 2;
      ctd2->CE[id].first++; 
      ctd2->CE[id].second += delE;
    }

    if (interFlag == ionize_2) {
      delE = IS.DIInterLoss(ch.IonList[Q2]);
      p2->iattrib[use_key] = k2.updateC(++C2);
      ctd2->CI[id].first++; 
      ctd2->CI[id].second += delE;
      partflag     = 2;
    }

    if (interFlag == recomb_2) {
      delE         = kEe1[id];	// See comment above for interFlag==6
      p2->iattrib[use_key] = k2.updateC(--C2);
      assert(C2 > 0);
      partflag     = 2;
      ctd2->RR[id].first++; 
      ctd2->RR[id].second += delE;
    }

    delEeV = delE;

    // Convert to super particle
    //
    if (partflag) delE *= NN;
    
    // Convert back to cgs
    //
    delE = delE * eV;
  }
  
  assert(delE >= 0.0);

  // Artifically prevent cooling
  //
  if (NO_COOL) delE = 0.0;

  // Pure scattering event
  //
  if (delE<=0.0) return ret;

  // Convert energy loss to system units
  //
  delE = delE/UserTreeDSMC::Eunit;
  
  // Assign interaction energy variables
  //
  double remE=0.0, totE=0.0, kEe=0.0, dof = 1.0;

  // Electrons from Particle 2 have interacted with atom/ion in Particle 1
  //
  if (partflag==1) {
    totE = kE;			// KE + internal
    if (use_Eint>=0) totE += p2->dattrib[use_Eint];
    remE = totE - dE;		// Energy floor
    kEe  = kEe2[id];		// Electron energy
    if (use_Eint>=0)
      dof = 1.0 + ne2;		// Total degrees of freedom

				// Energy diagnostics
    ctd1->eV_av[id] += kEe2[id];
    if (std::isnan(ctd2->eV_av[id])) {
      std::cout << "eV_N=" << ctd1->eV_N[id] << std::endl;
    }
    ctd1->eV_N[id]++;
    ctd1->eV_min[id] = std::min(ctd1->eV_min[id], kEe2[id]);
    ctd1->eV_max[id] = std::max(ctd2->eV_max[id], kEe2[id]);
    
    if (kEe2[id] > 10.2) { ctd1->eV_10[id]++;}
  }

  // Electrons from Particle 1 interacted with atom/ion in Particle 2
  //
  if (partflag==2) {
    totE = kE;			// KE + internal
    if (use_Eint>=0) totE += p1->dattrib[use_Eint];
    remE = totE - dE;		// Energy floor
    kEe  = kEe1[id];		// Electron energy
    if (use_Eint>=0)
      dof  = 1.0 + ne1;		// Total degrees of freedom

				// Energy diagnostics
    ctd2->eV_av[id] += kEe1[id];
    if (std::isnan(ctd2->eV_av[id])) {
      std::cout << "eV_N=" << ctd2->eV_N[id] << std::endl;
    }
    ctd2->eV_N[id]++;
    ctd2->eV_min[id] = std::min(ctd2->eV_min[id], kEe1[id]);
    ctd2->eV_max[id] = std::max(ctd2->eV_max[id], kEe1[id]);
    
    if (kEe1[id] > 10.2) { ctd2->eV_10[id]++; }
  }

  // Warn if energy lost is greater than total energy available to
  // lose
  //
  if (frost_warning && delE > remE)
      std::cout << "delE > KE!! (" << delE << " > " << totE
		<< "), Interaction type = " << interFlag 
		<< " kEe  = "  << kEe
		<< " delE = " << delEeV 
		<< std::endl;

  
  // Cooling rate diagnostic histogram
  //
  if (TSDIAG && delE>0.0) {
				// Histogram index
    int indx = (int)floor(log(remE/delE)/(log(2.0)*TSPOW) + 5);
				// Floor and ceiling
    if (indx<0 ) indx = 0;
    if (indx>10) indx = 10;
				// Add entry
    EoverT[id][indx] += Mt;
  }
  
  // Time step "cooling" diagnostic
  //
  if (use_delt>=0 && delE>0.0 && remE>0.0) {
    double dtE = delE/remE * spTau[id];
    double dt1 = p1->dattrib[use_delt];
    double dt2 = p2->dattrib[use_delt];
    p1->dattrib[use_delt] = std::max<double>(dt1, dtE);
    p2->dattrib[use_delt] = std::max<double>(dt2, dtE);
  }

  if (use_exes>=0) {
    // (-/+) value means under/overcooled: positive/negative increment
    // to delE NB: delE may be < 0 if too much energy was radiated
    // previously . . .
    //
    delE -= p1->dattrib[use_exes] + p2->dattrib[use_exes];
  }
  
  // Initial relative velocity 
  double vi        = (*cr);

  // Sufficient energy available for selected loss
  //
  if (remE > delE) {

    lostSoFar[id] += delE;
    decelT[id]    += delE;
    
    totE          -= delE;	// Remove the energy from the total
				// available

				// Energy per particle
    double kEm = totE / dof;
    
				// Get new relative velocity
    (*cr)          = sqrt( 2.0*kEm/Mu );

    ret            = 0;		// No error

    if (partflag==1) {
      ctd1->dv[id].first++; 
      ctd1->dv[id].second += 
	0.5*Mu*(vi - (*cr))*(vi - (*cr))/N1 * UserTreeDSMC::Eunit / eV;
    }
    
    if (partflag==2) {
      ctd2->dv[id].first++; 
      ctd2->dv[id].second += 
	0.5*Mu*(vi - (*cr))*(vi - (*cr))/N2 * UserTreeDSMC::Eunit / eV;
    }

				// Distribute electron energy to particles
    if (use_Eint>=0) {
      if (partflag==1) p2->dattrib[use_Eint] = ne2 * kEm;
      if (partflag==2) p1->dattrib[use_Eint] = ne1 * kEm;
    }
				// Zero out internal energy excess
    if (use_exes>=0)		// since excess is now used up
      p1->dattrib[use_exes] = p2->dattrib[use_exes] = 0.0;

  } else {
    //
    // Inconsistent: too much energy lost!
    //
    
    // Compute total available energy for both particles
    //
    totE = kE - dE;
    if (use_Eint>=0)
      totE += p1->dattrib[use_Eint] + p2->dattrib[use_Eint];

    // Try combined energy first . . . 
    //
    if (totE > delE) {

      lostSoFar[id] += delE;
      decelT[id]    += delE;
    
      totE          -= delE;	// Remove the energy from the total
				// available

				// Energy per particle
      double kEm = totE / (1.0 + ne1 + ne2);
    
				// Get new relative velocity
      (*cr)          = sqrt( 2.0*kEm/Mu );

      ret            = 0;	// No error

      if (partflag==1) {
	ctd1->dv[id].first++; 
	ctd1->dv[id].second += 
	  0.5*Mu*(vi - (*cr))*(vi - (*cr))/N1 * UserTreeDSMC::Eunit / eV;
      }
      
      if (partflag==2) {
	ctd2->dv[id].first++; 
	ctd2->dv[id].second +=
	  0.5*Mu*(vi - (*cr))*(vi - (*cr))/N2 * UserTreeDSMC::Eunit / eV;
      }

				// Remaining energy split between
				// internal degrees of freedom
      if (use_Eint>=0)
	p1->dattrib[use_Eint] = p2->dattrib[use_Eint] = kEm;
      
    } else {
				// All available energy will be lost
      lostSoFar[id] += totE;
      decolT[id]    += totE - delE;

      (*cr)         *= TolV;
      ret            = 1;	// Set error flag
    
				// Conservation of energy for internal
				// degrees of freedom
      dE             = 0.5*Mu*(*cr)*(*cr);

      if (partflag==1) {
	ctd1->dv[id].first++; 
	ctd1->dv[id].second +=
	  0.5*Mu*(vi - (*cr))*(vi - (*cr))/N1 * UserTreeDSMC::Eunit / eV;
      }
      
      if (partflag==2) {
	ctd2->dv[id].first++; 
	ctd2->dv[id].second +=
	  0.5*Mu*(vi - (*cr))*(vi - (*cr))/N2 * UserTreeDSMC::Eunit / eV;
      }

				// Remaining energy split set to zero
      if (use_Eint>=0)
	p1->dattrib[use_Eint] = p2->dattrib[use_Eint] = 0.0;

				// Reset internal energy excess
      if (use_exes>=0) {
	p1->dattrib[use_exes] = p1->mass*(totE - delE)/Mt;
	p2->dattrib[use_exes] = p2->mass*(totE - delE)/Mt;
      }
    }
  }
  
  return ret;
}


void CollideIon::debugDeltaE(double delE, unsigned short Z, unsigned short C,
			     double KE, double prob, int interFlag)
{
  std::cout << " *** Neg deltaE=" << std::setw(12) << delE << ", (Z, C)=(" 
	    << std::setw(2) << Z << ", " << std::setw(2) << C << "), E=" 
	    << std::setw(12) << KE  << ", prob=" << std::setw(12) << prob
	    << " :: " << labels[interFlag] << std::endl;
}

int CollideIon::inelasticTrace(pHOT *tree, Particle* p1, Particle* p2, 
			       double *cr, int id)
{
  // For particle map and weights
  //
  typedef std::map<speciesKey, double> keyW;
  typedef std::map<speciesKey, int   > keyI;

  int ret = 0;			// No error (flag)

  //-------------------------
  // VERBOSE DEBUG TEST
  //-------------------------
  // Set to false for production
  //          |
  //          v
  static bool DEBUG_F = false;
  //
  if (DEBUG_F) {
    static std::map<int, std::string> labels;
    if (labels.size()==0) {
      labels[geometric_1] = "geometric  [1]";
      labels[neut_elec_1] = "neutral el [1]";
      labels[ion_elec_1 ] = "charged el [1]";
      labels[free_free_1] = "free-free  [1]";
      labels[colexcite_1] = "col excite [1]";
      labels[ionize_1   ] = "ionization [1]";
      labels[recomb_1   ] = "recombine  [1]";
      labels[geometric_2] = "geometric  [2]";
      labels[neut_elec_2] = "neutral el [2]";
      labels[ion_elec_2 ] = "charged el [2]";
      labels[free_free_2] = "free-free  [2]";
      labels[colexcite_2] = "col excite [2]";
      labels[ionize_2   ] = "ionization [2]";
      labels[recomb_2   ] = "recombine  [2]";
    }
    //
    std::cout << std::setw( 8) << "index"
	      << std::setw( 4) << "Z1"
	      << std::setw( 4) << "C1"
	      << std::setw( 4) << "Z2"
	      << std::setw( 4) << "C2"
	      << std::setw( 8) << "flag"
	      << std::setw(14) << "cross"
	      << std::setw(18) << "type label"
	      << std::endl
	      << std::setw( 8) << "-----"
	      << std::setw( 4) << "---"
	      << std::setw( 4) << "---"
	      << std::setw( 4) << "---"
	      << std::setw( 4) << "---"
	      << std::setw( 8) << "-----"
	      << std::setw(14) << "---------"
	      << std::setw(18) << "---------------"
	      << std::endl;
    for (keyI::iterator sp=SpList.begin(); sp!=SpList.end(); sp++) {
      speciesKey k1(sp->first);
      for (keyI::iterator sp=SpList.begin(); sp!=SpList.end(); sp++) {
	speciesKey k2(sp->first);
	dKey dkey(k1, k2);
	for (size_t i = 0; i < sCrossMap[id][dkey].size(); i++) {
	  std::cout << std::setw( 8) << i
		    << std::setw( 4) << k1.first
		    << std::setw( 4) << k1.second
		    << std::setw( 4) << k2.first
		    << std::setw( 4) << k2.second
		    << std::setw( 8) << sInterMap[id][dkey][i]
		    << std::setw(14) << sCrossMap[id][dkey][i]
		    << std::setw(18) << labels[sInterMap[id][dkey][i]]
		    << std::endl;
	}
      }
    }
    std::cout << std::endl;
  }

  // Copy particle weights
  //
  keyW new1, new2;

  for (keyI::iterator sp=SpList.begin(); sp!=SpList.end(); sp++) {
    new1[sp->first] = p1->dattrib[sp->second];
    new2[sp->first] = p2->dattrib[sp->second];
  }

  // Cycle through all pairs of species
  //
  double delE = 0.0;

  keyCrossMap::iterator spBeg=sCrossMap[id].begin();
  keyCrossMap::iterator spEnd=sCrossMap[id].end();
  
  for (keyCrossMap::iterator sp=spBeg; sp!=spEnd; sp++) {

    // The interaction pair
    //
    dKey key = sp->first;

    // Number of interaction types in this map
    //
    size_t snum = sInterMap[id][key].size();

    // Normalize cross section list
    //
    double norm = std::accumulate(sp->second.begin(), sp->second.end(), 0.0);
    
    // Sanity check: compute elastic and inelastic fractions
    //  |
    //  v
    if (0) {
      double elastic = 0.0, inelastic = 0.0;
      for (size_t i=0; i<snum; i++) {
	if (sInterMap[id][key][i] % 100 < 3)
	  elastic   += sCrossMap[id][key][i]/norm;
	else
	  inelastic += sCrossMap[id][key][i]/norm;
      }
    }

    // The interaction pair
    //
    speciesKey k1 = sp->first.first;
    speciesKey k2 = sp->first.second;
    
    collTDPtr ctd1 = (*collD)[k1], ctd2 = (*collD)[k2];

    // Get weights
    //
    double w1 = p1->dattrib[SpList[k1]];
    double w2 = p2->dattrib[SpList[k2]];

    // Get weighted mass
    //
    double m1 = p1->mass * w1;
    double m2 = p2->mass * w2;

    // Species keys
    //
    unsigned short Z1 = k1.first, C1 = k1.second;
    unsigned short Z2 = k2.first, C2 = k2.second;
  
    // Number of atoms in each super particle of trace type
    //
    double N1 = m1*UserTreeDSMC::Munit/(atomic_weights[k1.first]*amu);
    double N2 = m2*UserTreeDSMC::Munit/(atomic_weights[k2.first]*amu);
  
    // Cycle through the interaction list
    //
    for (size_t isp=0; isp<snum; isp++) {
      
      int interFlag = sInterMap[id][key][isp];

      double prob   = sCrossMap[id][key][isp] * spProb[id];

      // Particle 1 weight
      //
      double W1     = prob * atomic_weights[k1.first] / w1;

      // Particle 2 weight
      //
      double W2     = prob * atomic_weights[k2.first] / w2;
      
      // Accumulate the total energy lost in inelastic processes
      //
      double tdelE = 0.0;

      // Accumulate the total energy lost for each particle (diagnostic)
      //
      double delE1 = 0.0, delE2 = 0.0;

      // Indicate loss in particle 1 or 2 (diagnostic)
      //
      bool p1Flag = false, p2Flag = false;

      //--------------------------------------------------
      // Ion keys
      //--------------------------------------------------

      lQ Q1(Z1, C1), Q2(Z2, C2);

      //-------------------------
      // Particle 1 interactions
      //-------------------------

      if (interFlag == free_free_1) {
	delE1 = IS.selectFFInteract(ch.IonList[Q1], kEe2[id]) * prob;
	if (delE1<0) debugDeltaE(delE1, Z1, C1, kEe2[id], prob, interFlag);
	tdelE               += delE1;
	ctd1->ff[id].first  += prob;
	ctd1->ff[id].second += delE1 * N1;
	p1Flag = true;
      }

      if (interFlag == colexcite_1) {
	delE1 = IS.selectCEInteract(ch.IonList[Q1], kCE1[id][key]) * prob;
	if (delE1<0) debugDeltaE(delE1, Z1, C1, 
				 kCE1[id][key].back().second, prob, interFlag);
	tdelE               += delE1;
	ctd1->CE[id].first  += prob;
	ctd1->CE[id].second += delE1 * N1;
	p1Flag = true;
      }

      if (interFlag == ionize_1) {
	delE1  = IS.DIInterLoss(ch.IonList[Q1]) * prob;
	if (delE1<0) debugDeltaE(delE1, Z1, C1, 0.0, prob, interFlag);
	tdelE += delE1;
	speciesKey kk(Z1, ++C1);
	if (W1 < (w1=new1[k1])) {
	  new1[kk] += W1;
	  new1[k1] -= W1;
	} else {
	  new1[kk] += w1;
	  new1[k1]  = 0.0;
	}
	assert(C1 <= Z1+1);
	ctd1->CI[id].first  += prob;
	ctd1->CI[id].second += delE1 * N1;
	p1Flag = true;
      }

      // KE carried by electron is subtracted from the thermal reservoir
      // but not the radiation of "binding".  The "binding" radiation
      // decreases the total energy of the gas but not the thermal
      // component.
      //
      if (interFlag == recomb_1) {
	delE1  = kEe2[id] * prob;
	if (delE1<0) debugDeltaE(delE1, Z1, C1, kEe2[id], prob, interFlag);
	tdelE += delE1;
	speciesKey kk(Z1, --C1);
	if (W1 < (w1=new1[k1])) {
	  new1[kk] += W1;
	  new1[k1] -= W1;
	} else {
	  new1[kk] += w1;
	  new1[k1]  = 0.0;
	}
	assert(C1 > 0);
	ctd1->RR[id].first  += prob;
	ctd1->RR[id].second += delE1 * N1;
	p1Flag = true;
      }
    
      //-------------------------
      // Particle 2 interactions
      //-------------------------
      
      if (interFlag == free_free_2) {
	delE2                = IS.selectFFInteract(ch.IonList[Q2], kEe1[id]) * prob;
	if (delE2<0) debugDeltaE(delE2, Z2, C2, kEe1[id], prob, interFlag);
	tdelE               += delE2;
	ctd2->ff[id].first  += prob;
	ctd2->ff[id].second += delE2 * N2;
	p2Flag = true;
      }

      if (interFlag == colexcite_2) {
	delE2 = IS.selectCEInteract(ch.IonList[Q2], kCE2[id][key]) * prob;
	if (delE2<0) debugDeltaE(delE2, Z2, C2, 
				 kCE2[id][key].back().second, prob, interFlag);
	tdelE               += delE2;
	ctd2->CE[id].first  += prob;
	ctd2->CE[id].second += delE2 * N2;
	p2Flag = true;
      }

      if (interFlag == ionize_2) {
	delE2   = IS.DIInterLoss(ch.IonList[Q2]) * prob;
	if (delE2<0) debugDeltaE(delE2, Z2, C2, 0.0, prob, interFlag);
	tdelE  += delE2;
	speciesKey kk(Z2, ++C2);
	if (W2 < (w2=new2[k2])) {
	  new2[kk] += W2;
	  new2[k2] -= W2;
	} else {
	  new2[kk] += w2;
	  new2[k2]  = 0.0;
	}
	assert(C2 <= Z2+1);
	ctd2->CI[id].first  += prob;
	ctd2->CI[id].second += delE2 * N2;
	p2Flag = true;
      }

      if (interFlag == recomb_2) {
	delE2  = kEe1[id] * prob; // See comment above for interFlag==6
	if (delE2<0) debugDeltaE(delE2, Z2, C2, kEe1[id], prob, interFlag);
	tdelE += delE2;
	speciesKey kk(Z2, --C2);
	if (W2 < (w2=new2[k2])) {
	  new2[kk] += W2;
	  new2[k2] -= W2;
	} else {
	  new2[kk] += w2;
	  new2[k2]  = 0.0;
	}
	assert(C2 > 0);
	ctd2->RR[id].first  += prob;
	ctd2->RR[id].second += delE2 * N2;
	p2Flag = true;
      }

      // Convert to super particle
      //
      if (interFlag/100 == 1) tdelE *= N1;
      if (interFlag/100 == 2) tdelE *= N2;
    
      // Energy diagnostics
      //
      if (p1Flag) {
	double kEe = kEe1[id];

	ctd1->eV_av[id] += kEe;
	if (std::isnan(ctd2->eV_av[id])) {
	  std::cout << "eV_N=" << ctd1->eV_N[id] << std::endl;
	}
	ctd1->eV_N[id]++;
	ctd1->eV_min[id] = std::min(ctd1->eV_min[id], kEe);
	ctd1->eV_max[id] = std::max(ctd2->eV_max[id], kEe);
	
	if (kEe > 10.2) { ctd1->eV_10[id]++;}

	ctd1->dv[id].first++; 
	ctd1->dv[id].second += delE1;
      }


      if (p2Flag) {
	double kEe = kEe2[id];

	ctd2->eV_av[id] += kEe;
	if (std::isnan(ctd2->eV_av[id])) {
	  std::cout << "eV_N=" << ctd2->eV_N[id] << std::endl;
	}
	ctd2->eV_N[id]++;
	ctd2->eV_min[id] = std::min(ctd2->eV_min[id], kEe);
	ctd2->eV_max[id] = std::max(ctd2->eV_max[id], kEe);
	
	if (kEe > 10.2) { ctd2->eV_10[id]++;}

	ctd1->dv[id].first++; 
	ctd1->dv[id].second += delE2;
      }
      
      // Convert back to cgs
      //
      delE += tdelE * eV;

      if (delE < 0.0) {
	std::cout << "Found delE=" << std::setw(14) << delE/UserTreeDSMC::Eunit 
		  << ", tdelE=" << std::setw(14) << tdelE
		  << ", delE1=" << std::setw(14) << delE1
		  << ", delE2=" << std::setw(14) << delE2
		  << ", w1="    << std::setw(14) << w1
		  << ", w2="    << std::setw(14) << w2
		  << ", N1="    << std::setw(14) << N1
		  << ", N2="    << std::setw(14) << N2
		  << ", (Z1, C1) = (" 
		  << std::setw(2) << Z1 << ", "
		  << std::setw(2) << C1 << ") "
		  << ", (Z2, C2) = (" 
		  << std::setw(2) << Z2 << ", "
		  << std::setw(2) << C2 << ") "
		  << std::endl;
      }
    }
  }
  
  if (delE < 0.0) {
    std::cout << "Found delE=" << delE/UserTreeDSMC::Eunit 
	      << " < 0.0" << std::endl;
    delE = 0.0;
  }

  // The total mass in system units
  //
  double m1 = p1->mass;
  double m2 = p2->mass;
  double Mt = m1 + m2;

  // Reduced mass in ballistic collision (system units)
  //
  double Mu = m1 * m2 / Mt;

  // Center of mass energy in the ballistic collision (system units)
  //
  double kE  = 0.5*Mu*(*cr)*(*cr);

  // For tracking energy conservation (system units)
  //
  double dE   = kE*TolV*TolV;
    
  // Artifically prevent cooling
  //
  if (NO_COOL) delE = 0.0;

  // Convert energy loss to system units
  //
  delE = delE/UserTreeDSMC::Eunit;
  
  // Assign interaction energy variables
  //
  double remE, totE, kEe = kEe1[id];

  // Diagnostic accumulation
  //
  totE = kE;			// KE
  remE = totE - dE;		// Energy floor
  
  // Warn if energy lost is greater than total energy available to
  // lose
  //
  if (frost_warning && delE > remE)
    std::cout << "delE > KE!! (" << delE << " > " << totE
	      << "), kEe  = "  << kEe
	      << " delE = " << delE/(eV*Mu*UserTreeDSMC::Munit*amu)
	      << std::endl;
  
  // Cooling rate diagnostic histogram
  //
  if (TSDIAG && delE>0.0) {
				// Histogram index
    int indx = (int)floor(log(remE/delE)/(log(2.0)*TSPOW) + 5);
				// Floor and ceiling
    if (indx<0 ) indx = 0;
    if (indx>10) indx = 10;
				// Add entry
    EoverT[id][indx] += Mt;
  }
  
  // Time step "cooling" diagnostic
  //
  if (use_delt>=0 && delE>0.0 && remE>0.0) {
    double dtE = delE/remE * spTau[id];
    double dt1 = p1->dattrib[use_delt];
    double dt2 = p2->dattrib[use_delt];
    p1->dattrib[use_delt] = std::max<double>(dt1, dtE);
    p2->dattrib[use_delt] = std::max<double>(dt2, dtE);
  }

  // Inelastic
  //
  if (delE>0.0) {

    // Sufficient energy available for selected loss
    //
    if (remE > delE) {
      
      lostSoFar[id] += delE;
      decelT[id]    += delE;
      
      totE          -= delE;	// Remove the energy from the total
				// available
      
				// Energy per particle
    
				// Get new relative velocity
      (*cr)          = sqrt( 2.0*totE/Mu );

      ret            = 0;		// No error

    } else {
      //
      // Inconsistent: too much energy lost!
      //
      
      lostSoFar[id] += totE;
      decolT[id]    += totE - delE;

      (*cr)         *= TolV;
      ret            = 1;	// Set error flag
      
				// Conservation of energy for internal
				// degrees of freedom
      dE             = 0.5*Mu*(*cr)*(*cr);

    }
  }
  
  // Replace particle weights
  //
  for (keyI::iterator sp=SpList.begin(); sp!=SpList.end(); sp++) {
    p1->dattrib[sp->second] = new1[sp->first];
    p2->dattrib[sp->second] = new2[sp->first];
  }

  // Velocity update
  //
  std::vector<double> vrel_i(3), vrel_f(3), vcm(3);

  // Pre-collision relative velocity
  // 
  for(unsigned k=0; k<3; k++)
    vrel_i[k] = p1->vel[k] - p2->vel[k];

  for(unsigned k=0; k<3; k++)
    vcm[k] = (m1*p1->vel[k] + m2*p2->vel[k]) / Mt;
	    
  double cos_th = 1.0 - 2.0*(*unit)();       // Cosine and sine of
  double sin_th = sqrt(1.0 - cos_th*cos_th); // collision angle theta
  double phi    = 2.0*M_PI*(*unit)();	     // Collision angle phi
  
  vrel_f[0] = (*cr)*cos_th;	     // Compute post-collision
  vrel_f[1] = (*cr)*sin_th*cos(phi); // relative velocity
  vrel_f[2] = (*cr)*sin_th*sin(phi);
  
  // Update post-collision velocities
  // 
  for (unsigned k=0; k<3; k++) {
    p1->vel[k] = vcm[k] + m2/Mt*vrel_f[k];
    p2->vel[k] = vcm[k] - m1/Mt*vrel_f[k];
  }

  return ret;
}


double CollideIon::Etotal()
{ 
  double ret = totalSoFar;
  totalSoFar = 0.0;
  return ret; 
}

double CollideIon::Mtotal()
{ 
  double ret = massSoFar;
  massSoFar = 0.0;
  return ret; 
}

void CollideIon::Elost(double* collide, double* epsm)
{ 
  double ret1=0.0;
  for (int n=0; n<nthrds; n++) {
    ret1 += lostSoFar[n];
    lostSoFar[n] = 0.0; 
  }
  *collide = ret1;
}


void * CollideIon::timestep_thread(void * arg)
{
  pHOT* tree = (pHOT* )((tstep_pass_arguments*)arg)->tree;
  int id     = (int)((tstep_pass_arguments*)arg)->id;
  
  thread_timing_beg(id);
  
  // Loop over cells, cell time-of-flight time
  // for each particle
  
  pCell *c;
  Particle *p;
  double DT, mscale;
  // double L;
  
  // Loop over cells, processing collisions in each cell
  //
  for (unsigned j=0; j<cellist[id].size(); j++ ) {
    
    // Number of particles in this cell
    //
    c = cellist[id][j];
    
    double volc = c->Volume();
    
    sKeyDmap            densM, lambdaM, crossM;
    sKeyUmap::iterator  it1, it2;
    sKey2Dmap           crossIJ;
    
    crossIJ = totalScatteringCrossSections(0, c, id);
    
    for (it1=c->count.begin(); it1!=c->count.end(); it1++) {
      speciesKey i1 = it1->first;
      densM[i1] = c->Mass(i1)/volc;
    }
    
    double meanDens=0.0, meanLambda=0.0;
    
    for (it1=c->count.begin(); it1!=c->count.end(); it1++) {
      speciesKey i1 = it1->first;
      crossM [i1] = 0.0;
      for (it2=c->count.begin(); it2!=c->count.end(); it2++) {

	speciesKey i2 = it2->first;
	double      N = UserTreeDSMC::Munit/(amu*atomic_weights[i2.first]);
	
	if (i2>=i1) {
	  crossM[i1] += N*densM[i2] * crossIJ[i1][i2];
	} else
	  crossM[i1] += N*densM[i2] * crossIJ[i2][i1];
      }
      
      lambdaM[i1] = 1.0/crossM[i1];
      meanDens   += densM[i1] ;
      meanLambda += densM[i1] * lambdaM[i1];
    }
    
    // This is the number density-weighted
    meanLambda /= meanDens;
    
    for (vector<unsigned long>::iterator 
	   i=c->bods.begin(); i!=c->bods.end(); i++) {
      // Current particle
      p = tree->Body(*i);
      // Compute time of flight criterion
      DT = 1.0e40;
      mscale = 1.0e40;
      for (unsigned k=0; k<3; k++) {
	DT = min<double>
	  (meanLambda/(fabs(p->vel[k])+1.0e-40), DT);
      }
      mscale = min<double>(meanLambda, mscale);

      // Size scale for multistep timestep calc.
      p->scale = mscale;

      // Compute cooling criterion timestep
      if (use_delt>=0) {
	double v = p->dattrib[use_delt];
	if (v>0.0) DT = min<double>(DT, v);
      }
      
      p->dtreq = DT;
    }
  }
  
  thread_timing_end(id);
  
  return (NULL);
}

void CollideIon::finalize_cell(pHOT* tree, pCell* cell, double kedsp, int id)
{
  //
  // Cross-section debugging [END]
  //
  if (CROSS_DBG && id==0) {
    if (nextTime_dbg <= tnow && nCnt_dbg < nCel_dbg)
      {
	speciesKey i;

	if (aType==Direct)
	  i = cell->count.begin()->first;
	else
	  i = SpList.begin()->first;

	cross2_dbg.push_back(csections[id][i][i]);
	nCnt_dbg++;
	if (nCnt_dbg == nCel_dbg) write_cross_debug();
      }
  }
  // Done

}

// Help class that maintains database of diagnostics
//
collDiag::collDiag(CollideIon* caller) : p(caller)
{
  // Initialize the map
  //
  if (p->ZList.size()) {
    for (CollideIon::ZLtype::iterator 
	   n=p->ZList.begin(); n!=p->ZList.end(); n++) {

      unsigned short Z = *n;

      for (unsigned short C=1; C<Z+2; C++) {
	speciesKey k(Z, C);
	(*this)[k] = collTDPtr(new CollisionTypeDiag());
      }
    }
  } else if (p->SpList.size()) {
    typedef std::map<speciesKey, int>::iterator spMapItr;
    for (spMapItr n=p->SpList.begin(); n!=p->SpList.end(); n++) {
      (*this)[n->first] = collTDPtr(new CollisionTypeDiag());
    }
  } else {			// Sanity check
    if (myid==0) {
      std::cerr << "collDiag:CollDiag: species list or map is "
		<< "not initialized" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 57);
  }

  // Initialize the output file
  //
  initialize();
}

// Gather statistics from all processes
//
void collDiag::gather()
{
  for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
    collTDPtr ctd = it->second;
    ctd->sumUp();

    double z;
    if (myid==0) {
      MPI_Reduce(&(z=ctd->ff_s.first),  &ctd->ff_s.first,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 1
      MPI_Reduce(&(z=ctd->ff_s.second), &ctd->ff_s.second, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 2
      MPI_Reduce(&(z=ctd->CE_s.first),  &ctd->CE_s.first,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 3
      MPI_Reduce(&(z=ctd->CE_s.second), &ctd->CE_s.second, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 4
      MPI_Reduce(&(z=ctd->CI_s.first),  &ctd->CI_s.first,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 5
      MPI_Reduce(&(z=ctd->CI_s.second), &ctd->CI_s.second, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 6
      MPI_Reduce(&(z=ctd->RR_s.first),  &ctd->RR_s.first,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 7
      MPI_Reduce(&(z=ctd->RR_s.second), &ctd->RR_s.second, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 8
      MPI_Reduce(&(z=ctd->dv_s.first),  &ctd->dv_s.first,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 9
      MPI_Reduce(&(z=ctd->dv_s.second), &ctd->dv_s.second, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 10
      MPI_Reduce(&(z=ctd->eV_av_s),     &ctd->eV_av_s,     1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 11
      MPI_Reduce(&(z=ctd->eV_N_s),      &ctd->eV_N_s,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 12
      MPI_Reduce(&(z=ctd->eV_min_s),    &ctd->eV_min_s,    1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); // 13
      MPI_Reduce(&(z=ctd->eV_max_s),    &ctd->eV_max_s,    1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); // 14
      MPI_Reduce(&(z=ctd->eV_10_s),     &ctd->eV_10_s,     1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 15
    } else {
      MPI_Reduce(&ctd->ff_s.first,      &z,                1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 1 
      MPI_Reduce(&ctd->ff_s.second, 	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 2 
      MPI_Reduce(&ctd->CE_s.first,  	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 3 
      MPI_Reduce(&ctd->CE_s.second, 	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 4 
      MPI_Reduce(&ctd->CI_s.first,  	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 5 
      MPI_Reduce(&ctd->CI_s.second, 	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 6 
      MPI_Reduce(&ctd->RR_s.first,  	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 7 
      MPI_Reduce(&ctd->RR_s.second, 	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 8 
      MPI_Reduce(&ctd->dv_s.first,    	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 9 
      MPI_Reduce(&ctd->dv_s.second,   	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 10
      MPI_Reduce(&ctd->eV_av_s,       	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 11
      MPI_Reduce(&ctd->eV_N_s,        	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 12
      MPI_Reduce(&ctd->eV_min_s,      	&z, 		   1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD); // 13
      MPI_Reduce(&ctd->eV_max_s,      	&z, 		   1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); // 14
      MPI_Reduce(&ctd->eV_10_s,       	&z, 		   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 15
    }
  }
}

// Zero out the counters
//
void collDiag::reset()
{
  for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++)  
    it->second->reset();
}

void collDiag::initialize() 
{
  if (myid) return;

  {
    // Generate the file name
    std::ostringstream sout;
    sout << outdir << runtag << ".ION_coll_debug";
    coll_file_debug = sout.str();

    // Check for existence
    std::ifstream in(coll_file_debug.c_str());

    if (in.fail()) {
      // Write a new file
      std::ofstream out(coll_file_debug.c_str());
      if (out) {

	out << "# Variable      key                      " << std::endl
	    << "# ------------  -------------------------" << std::endl
	    << "# N(ff)         number of free-free      " << std::endl
	    << "# E(ff)         cum energy in free-free  " << std::endl
	    << "# N(ce)         number of collions       " << std::endl
	    << "# E(ce)         cum energy in collisions " << std::endl
	    << "# N(ci)         number of ionizations    " << std::endl
	    << "# E(ci)         cum energy in ionizations" << std::endl
	    << "# N(rr)         number of rad recombs    " << std::endl
	    << "# E(rr)         energy in rad recombs    " << std::endl
	    << "# d(KE)         mean energy change       " << std::endl
	    << "# Etotl         total cum energy         " << std::endl
	    << "#"                                         << std::endl;
	
				// Species labels
	out << "#" << std::setw(11) << "Species==>" << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  ostringstream sout, sout2;
	  sout  << "(" << it->first.first << ", " << it->first.second << ")";
	  size_t w = 9*12, l = sout.str().size();
	  sout2 << std::setw((w-l)/2) << ' ' << sout.str();
	  out   << std::setw(w) << sout2.str() << " | ";
	}
	out << std::setw(12) << ' ' << " |" << std::endl;

				// Header line
	out << std::setfill('-') << std::right;
	out << "#" << std::setw(11) << '+' << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  for (int i=0; i<9; i++) out << std::setw(12) << '+';
	  out << " | ";
	}
	out << std::setw(12) << '-' << " |"
	    << std::setfill(' ') << std::endl;

				// Column labels
	out << "#" << std::setw(11) << "Time" << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  out << std::setw(12) << "N(ff) |"
	      << std::setw(12) << "E(ff) |"
	      << std::setw(12) << "N(ce) |"
	      << std::setw(12) << "E(ce) |"
	      << std::setw(12) << "N(ci) |"
	      << std::setw(12) << "E(ci) |"
	      << std::setw(12) << "N(rr) |"
	      << std::setw(12) << "E(rr) |"
	      << std::setw(12) << "d(KE) |"
	      << " | ";
	}
	out << std::setw(12) << "Etotl" << std::endl;

				// Column indices
	unsigned cnt = 0;
	out << "#" << std::setw(11) << std::right << ++cnt << " | ";
	std::ostringstream ts;
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  for (size_t l=0; l<9; l++) {
	    ts.str("");
	    ts << "[" << ++cnt << "] |";
	    out << std::setw(12) << std::right << ts.str();
	  }
	  out << " | ";
	}
	ts.str("");
	ts << ++cnt << " |";
	out << std::setw(12) << std::right << ts.str() << std::endl;

				// Header line
	out << std::setfill('-') << std::right;
	out << "#" << std::setw(11) << '+' << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  for (int i=0; i<9; i++) out << std::setw(12) << '+';
	  out << " | ";
	}
	out << std::setw(12) << '-' << " |"
	    << std::setfill(' ') << std::endl;
      }
    }
    in.close();
  }

  {
    // Generate the file name
    std::ostringstream sout;
    sout << outdir << runtag << ".ION_energy_debug";
    energy_file_debug = sout.str();
    
    // Check for existence
    std::ifstream in(energy_file_debug.c_str());
    
    if (in.fail()) {
      // Write a new file
      std::ofstream out(energy_file_debug.c_str());
      if (out) {
	
	out << "# Variable      key                      " << std::endl
	    << "# ------------  -------------------------" << std::endl
	    << "# avg           mean collision energy    " << std::endl
	    << "# num           number of collisions     " << std::endl
	    << "# min           minimum collison energy  " << std::endl
	    << "# max           maximum collision energy " << std::endl
	    << "# over10        number > 10.2 eV         " << std::endl
	    << "#"                                         << std::endl;

				// Species labels
	out << "#" << std::setw(11) << "Species==>" << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  ostringstream sout, sout2;
	  sout  << "(" << it->first.first << ", " << it->first.second << ")";
	  size_t w = 5*12, l = sout.str().size();
	  sout2 << std::setw((w-l)/2) << ' ' << sout.str();
	  out   << std::setw(w) << sout2.str() << " | ";
	}
	out << std::endl;
	
				// Header line
	out << std::setfill('-') << std::right;
	out << "#" << std::setw(11) << '+' << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  for (int i=0; i<5; i++) out << std::setw(12) << '+';
	  out << " | ";
	}
	out << std::setfill(' ') << std::endl;

				// Column labels
	out << "#" << std::setw(11) << "Time" << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  out << std::setw(12) << "avg |"
	      << std::setw(12) << "num |"
	      << std::setw(12) << "min |"
	      << std::setw(12) << "max |"
	      << std::setw(12) << "over10 |"
	      << " | ";
	}
	out << std::endl;

				// Header line
	out << std::setfill('-') << std::right;
	out << "#" << std::setw(11) << '+' << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  for (int i=0; i<5; i++) out << std::setw(12) << '+';
	  out << " | ";
	}
	out << std::setfill(' ') << std::endl;
      }
    }
    in.close();
  }

}

void collDiag::print() 
{
  if (myid) return;

  {
    std::ofstream out(coll_file_debug.c_str(), ios::out | ios::app);
    out << std::scientific << std::setprecision(3);
    if (out) {
      double Etot = 0.0;
      double cvrt = eV/UserTreeDSMC::Eunit;
      out << std::setw(12) << tnow << " | ";
      for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	collTDPtr ctd = it->second;
	out << std::setw(12) << ctd->ff_s.first
	    << std::setw(12) << ctd->ff_s.second * cvrt
	    << std::setw(12) << ctd->CE_s.first
	    << std::setw(12) << ctd->CE_s.second * cvrt
	    << std::setw(12) << ctd->CI_s.first
	    << std::setw(12) << ctd->CI_s.second * cvrt
	    << std::setw(12) << ctd->RR_s.first
	    << std::setw(12) << ctd->RR_s.second * cvrt;
	if (ctd->dv_s.first>0.0)
	  out << std::setw(12) << ctd->dv_s.second/ctd->dv_s.first << " | ";
	else
	  out << std::setw(12) << 0.0 << " | ";
	Etot += 
	  ctd->ff_s.second + ctd->CE_s.second +
	  ctd->CI_s.second + ctd->RR_s.second ;
      }
      out << std::setw(12) << Etot * cvrt << " |" << std::endl;
    }
  }

  {
    std::ofstream out(energy_file_debug.c_str(), ios::out | ios::app);
    if (out) {
      out << std::setw(12) << tnow << " | ";
      for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	collTDPtr ctd = it->second;
	out << std::setw(12) << (ctd->eV_N_s ? ctd->eV_av_s/ctd->eV_N_s : 0.0)
	    << std::setw(12) << ctd->eV_N_s
	    << std::setw(12) << ctd->eV_min_s
	    << std::setw(12) << ctd->eV_max_s
	    << std::setw(12) << ctd->eV_10_s
	    << " | ";
      }
      out << std::endl;
    }
  }

}

void CollideIon::parseSpecies(const std::string& map)
{
  std::ifstream in(map.c_str());
  if (in.bad()) 
    {
      if (myid==0) {
	std::cerr << "CollideIon::parseSpecies: definition file <"
		  << map << "> could not be opened . . . quitting"
		  << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 55);
    }


  // Ok, read the first line to get the implementation type

  const int nline = 2048;
  char line[nline];

  in.getline(line, nline);
  std::string type(line);

  if (type.compare("direct")==0) {
    
    aType = Direct;

    if (use_key<0) {
      if (myid==0) {
	std::cerr << "UserTreeDSMC: species key position is not defined in "
		  << "Component" << std::endl;
      }
      MPI_Abort(MPI_COMM_WORLD, 36);
    }

    int Z;
    while (1) {
      in.getline(line, nline);
      if (in.good()) {
	std::istringstream sz(line);
	sz >> Z;		// Add to the element list
	if (!sz.bad()) ZList.insert(Z);
      } else {
	break;
      }
    }

  } else if (type.compare("trace")==0) {

    aType = Trace;
    
    speciesKey key;
    int pos;
    while (1) {
      in.getline(line, nline);
      if (in.good()) {
	std::istringstream sz(line);
	sz >> key.first;
	sz >> key.second;
	sz >> pos;
				// Add to the species list
	if (!sz.bad()) {
	  SpList[key] = pos;
	  ZList.insert(key.first);
	}
      } else {
	break;
      }
    }

  } else {
    if (myid==0) {
      std::cerr << "CollideIon::parseSpecies: implementation type <"
		<< type << "> is not recognized . . . quitting"
		<< std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 56);
  }

}

sKey2Umap CollideIon::generateSelection
(pCell* c, sKeyDmap* Fn, double crm, double tau, int id,
 double& meanLambda, double& meanCollP, double& totalNsel)
{
  if (aType == Direct)
    return generateSelectionDirect(c, Fn, crm, tau, id, 
				   meanLambda, meanCollP, totalNsel);
  else
    return generateSelectionTrace(c, Fn, crm, tau, id, 
				  meanLambda, meanCollP, totalNsel);
}

sKey2Umap CollideIon::generateSelectionDirect
(pCell* c, sKeyDmap* Fn, double crm, double tau, int id,
 double& meanLambda, double& meanCollP, double& totalNsel)
{
  sKeyDmap            densM, collPM, lambdaM, crossM;
  sKey2Dmap           selcM;
  sKey2Umap           nselM;
  sKeyUmap::iterator  it1, it2;
    
  // Volume in the cell
  //
  double volc = c->Volume();
  
  //
  // Cross-section debugging [BEGIN]
  //
  
  if (CROSS_DBG && id==0) {
    if (nextTime_dbg <= tnow && nCnt_dbg < nCel_dbg) {
      speciesKey i = c->count.begin()->first;
      cross1_dbg.push_back(csections[id][i][i]);
    }
  }
  // Done
  
  for (it1=c->count.begin(); it1!=c->count.end(); it1++) {
    speciesKey i1 = it1->first;
    densM[i1] = c->Mass(i1)/volc;
  }
    
  double meanDens = 0.0;
  meanLambda      = 0.0;
  meanCollP       = 0.0;
    
  for (it1=c->count.begin(); it1!=c->count.end(); it1++) {

    speciesKey i1 = it1->first;
    crossM [i1]   = 0.0;

    for (it2=c->count.begin(); it2!=c->count.end(); it2++) {
      speciesKey i2 = it2->first;

      if (i2>=i1) {
	crossM[i1] += (*Fn)[i2]*densM[i2]*csections[id][i1][i2];
      } else
	crossM[i1] += (*Fn)[i2]*densM[i2]*csections[id][i2][i1];
      
      if (csections[id][i1][i2] <= 0.0 || isnan(csections[id][i1][i2])) {
	cout << "INVALID CROSS SECTION! :: " << csections[id][i1][i2]
	     << " #1 = (" << i1.first << ", " << i1.second << ")"
	     << " #2 = (" << i2.first << ", " << i2.second << ")";
	csections[id][i1][i2] = 0.0; // Zero out
      }
	    
      if (csections[id][i2][i1] <= 0.0 || isnan(csections[id][i2][i1])) {
	cout << "INVALID CROSS SECTION! :: " << csections[id][i2][i1]
	     << " #1 = (" << i2.first << ", " << i2.second << ")"
	     << " #2 = (" << i1.first << ", " << i1.second << ")";
	csections[id][i2][i1] = 0.0; // Zero out
      }
	
    }
      
    if (it1->second>0 && (crossM[i1] == 0 || isnan(crossM[i1]))) {
      cout << "INVALID CROSS SECTION! ::"
	   << " crossM = " << crossM[i1] 
	   << " densM = "  <<  densM[i1] 
	   << " Fn = "     <<  (*Fn)[i1] << endl;
    }
    
    lambdaM[i1] = 1.0/crossM[i1];
    collPM [i1] = crossM[i1] * crm * tau;
    
    meanDens   += densM[i1] ;
    meanCollP  += densM[i1] * collPM[i1];
    meanLambda += densM[i1] * lambdaM[i1];
  }
    
  // This is the number density-weighted
  // MFP (used for diagnostics only)
  meanLambda /= meanDens;

  // Number-density weighted collision
  // probability (used for diagnostics
  // only)
  meanCollP  /= meanDens;
    
  // This is the per-species N_{coll}
  //
  totalNsel = 0.0;
  for (it1=c->count.begin(); it1!=c->count.end(); it1++) {
    speciesKey i1 = it1->first;
    
    for (it2=it1; it2!=c->count.end(); it2++) {
      speciesKey i2 = it2->first;
      
      // Probability of an interaction of between particles of type 1
      // and 2 for a given particle of type 2
      double Prob = (*Fn)[i2] * densM[i2] * csections[id][i1][i2] * crm * tau;
      
      if (i1==i2)
	selcM[i1][i2] = 0.5 * (it1->second-1) *  Prob;
      else
	selcM[i1][i2] = it1->second * Prob;
      
      nselM[i1][i2] = static_cast<unsigned>(floor(selcM[i1][i2]+0.5));
      totalNsel += nselM[i1][i2];
    }
  }
  
  return nselM;
}

sKey2Umap CollideIon::generateSelectionTrace
(pCell* c, sKeyDmap* Fn, double crm, double tau, int id,
 double& meanLambda, double& meanCollP, double& totalNsel)
{
  speciesKey          key(defaultKey);
    
  // Mass density in the cell
  //
  double dens = c->Mass() / c->Volume();
  
  // Number of bodies in this cell
  //
  unsigned num = static_cast<unsigned>(c->bods.size());

  //
  // Cross-section debugging [BEGIN]
  //
  
  if (CROSS_DBG && id==0) {
    if (nextTime_dbg <= tnow && nCnt_dbg < nCel_dbg) {
      speciesKey i = c->count.begin()->first;
      cross1_dbg.push_back(csections[id][i][i]);
    }
  }
  // Done
  
  // Sanity check
  if (isnan(csections[id][key][key]) or csections[id][key][key] < 0.0) {
    cout << "INVALID CROSS SECTION! :: " << csections[id][key][key] << std::endl;
    csections[id][key][key] = 0.0; // Zero out
  }
    
  // Cache relative velocity
  //
  spCrm[id] = crm;

  // Cross sections are already number weighted at this point
  //
  double crossM = (*Fn)[key] * dens * csections[id][key][key];
  double collPM = crossM * crm * tau;

  // Cache probability of an interaction of between the particles pair
  // for use in inelasticTrace
  //
  double rateF = 0.5 * (*Fn)[key] * dens * crm * tau;
  //             ^
  //             |
  // Cross sections are double counted

  spProb[id] = rateF * 1e-14 / (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);


  // Cache time step for estimating "over" cooling timestep is use_delt>=0
  //
  spTau[id]  = tau;

  // For Collide diagnostics
  //
  meanLambda = 1.0/crossM;
  meanCollP  = collPM;
    
  double Prob = rateF * csections[id][key][key];
  double selcM = 0.5 * (num-1) * Prob;
  //             ^
  //             |
  // Pairs are double counted

  sKey2Umap nselM;
  nselM[key][key] = static_cast<unsigned>(floor(selcM+0.5));
  totalNsel = selcM;
  
  return nselM;
}


void CollideIon::write_cross_debug()
{
  std::ofstream out(cross_debug.c_str(), ios::out | ios::app);
  for (int i=0; i<nCel_dbg; i++) {
    double diff = cross2_dbg[i] - cross1_dbg[i];
    if (cross1_dbg[i]>0.0) diff /= cross1_dbg[i];
    out << std::setw( 8) << i+1
	<< std::setw(18) << tnow
	<< std::setw(18) << cross1_dbg[i]
	<< std::setw(18) << cross2_dbg[i]
	<< std::setw(18) << diff
	<< std::endl;
  }
  nextTime_dbg += delTime_dbg;
  nCnt_dbg = 0;
  cross1_dbg.erase(cross1_dbg.begin(), cross1_dbg.end());
  cross2_dbg.erase(cross2_dbg.begin(), cross2_dbg.end());
}



void CollideIon::gatherSpecies()
{
  if (aType==Direct) return;

  // Trace version follows

  const double Tfac = 2.0*UserTreeDSMC::Eunit/3.0 * amu /
    UserTreeDSMC::Munit/boltz;

  // Clean the maps
  //
  double mass = 0.0;

  tempM = 0.0;
  specM.erase(specM.begin(), specM.end());

  // Interate through all cells
  //
  pHOT_iterator itree(*c0->Tree());

  while (itree.nextCell()) {

    pCell *cell = itree.Cell();

    // Compute the temerature
    //
    double KEtot, KEdsp;
    cell->sample->KE(defaultKey, KEtot, KEdsp);
    double T = KEdsp* Tfac;

    // Iterate through all bodies in this cell
    //
    vector<unsigned long>::iterator j = cell->bods.begin();
    while (j != cell->bods.end()) {
      Particle* p = cell->Body(j++);
      for (spItr it=SpList.begin(); it!=SpList.end(); it++) {

	speciesKey k = it->first;
	int     indx = it->second;
      
	if (specM.find(k) == specM.end()) specM[k] = 0.0;
	specM[k] += p->mass * p->dattrib[indx];
      } 
    }
    mass  += cell->Mass();
    tempM += cell->Mass() * T;
  }

  // Send the local map to other nodes
  //
  int sizm;
  spDItr it;
  speciesKey key;
  double val1, val2;

  for (int i=1; i<numprocs; i++) {
    if (i == myid) {
      sizm = specM.size();
				// Local map size
      MPI_Send(&sizm,  1, MPI_INT,    0, 330, MPI_COMM_WORLD);
				// Mass
      MPI_Send(&mass,  1, MPI_DOUBLE, 0, 331, MPI_COMM_WORLD);
				// Temp
      MPI_Send(&tempM, 1, MPI_DOUBLE, 0, 332, MPI_COMM_WORLD);
				// Send local map
      for (it=specM.begin(); it != specM.end(); it++) {
	key  = it->first;
	MPI_Send(&key.first,  1, MPI_UNSIGNED_SHORT, 0, 333, MPI_COMM_WORLD);
	MPI_Send(&key.second, 1, MPI_UNSIGNED_SHORT, 0, 334, MPI_COMM_WORLD);
	MPI_Send(&it->second, 1, MPI_DOUBLE,         0, 335, MPI_COMM_WORLD);
      }

    }
				// Root receives from Node i
    if (0 == myid) {
      MPI_Recv(&sizm, 1, MPI_INT,    i, 330, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&val1, 1, MPI_DOUBLE, i, 331, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&val2, 1, MPI_DOUBLE, i, 332, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      mass  += val1;
      tempM += val2;
      for (int j=0; j<sizm; j++) {
	MPI_Recv(&key.first,  1, MPI_UNSIGNED_SHORT, i, 333, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	MPI_Recv(&key.second, 1, MPI_UNSIGNED_SHORT, i, 334, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	MPI_Recv(&val1,       1, MPI_DOUBLE,         i, 335, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
				// Update root's map
				// 
	if (specM.find(key) == specM.end()) specM[key] = 0.0;
	specM[key] += val1;
      }
    }
  }

				// At this point, root's map is global
				// and remaning nodes have local maps
  if (mass>0.0) {
    for (spDItr it=specM.begin(); it != specM.end(); it++) {
      it->second /= mass;
    }
    tempM /= mass;
  }
}
  

// Print out species counts
void CollideIon::printSpecies(std::map<speciesKey, unsigned long>& spec)
{
  if (myid) return;

  if (aType == Direct) {	// Call the generic printSpecies member
    Collide::printSpecies(spec);
  } else {			// Call the trace fraction version
    printSpeciesTrace();
  }

}

const std::string clabl(unsigned c)
{
  std::ostringstream sout;
  sout << "[" << c << "]  ";
  return sout.str();
}

void CollideIon::printSpeciesTrace()
{
  std::ofstream dout;

  // Generate the file name
  //
  if (species_file_debug.size()==0) {
    std::ostringstream sout;
    sout << outdir << runtag << ".species";
    species_file_debug = sout.str();

    // Open the file for the first time
    //
    dout.open(species_file_debug.c_str());

    // Print the header
    //
    dout << "# " 
	 << std::setw(12) << std::right << "Time  "
	 << std::setw(12) << std::right << "Temp  ";
    for (spDItr it=specM.begin(); it != specM.end(); it++) {
      std::ostringstream sout;
      sout << "(" << it->first.first << "," << it->first.second << ") ";
      dout << std::setw(12) << right << sout.str();
    }
    dout << std::endl;

    unsigned cnt = 0;
    dout << "# " 
	 << std::setw(12) << std::right << clabl(++cnt);
    dout << std::setw(12) << std::right << clabl(++cnt);
    for (spDItr it=specM.begin(); it != specM.end(); it++)
      dout << std::setw(12) << right << clabl(++cnt);
    dout << std::endl;

    dout << "# " 
	 << std::setw(12) << std::right << "--------"
	 << std::setw(12) << std::right << "--------";
    for (spDItr it=specM.begin(); it != specM.end(); it++)
      dout << std::setw(12) << std::right << "--------";
    dout << std::endl;

  } else {			// Open for append
				//
    dout.open(species_file_debug.c_str(), ios::out | ios::app);
  }

  dout << std::setprecision(5);
  dout << "  " 
       << std::setw(12) << std::right << tnow
       << std::setw(12) << std::right << tempM;
  for (spDItr it=specM.begin(); it != specM.end(); it++)
    dout << std::setw(12) << std::right << it->second;
  dout << std::endl;
}


// Compute the mean molecular weight in atomic mass units
double CollideIon::molWeight(Component *c0)
{
  // Only do this once and cache the result
  if (mol_weight<0)  {
    
    mol_weight = 1.0;		// Default; i.e. for trace algorithm

    if (aType==Direct) {

      mol_weight = 0.0;

      size_t j = 0;
      std::map<unsigned short, int> zval;
      for (ZLtype::iterator i=ZList.begin(); i!=ZList.end(); i++) 
	zval[*i] = j++;
      
      std::map<unsigned short, int>::iterator iz, izEnd(zval.end());
      
      std::vector<double> speciesM(ZList.size(), 0.0);
      PartMapItr p = c0->Particles().begin(), pEnd = c0->Particles().end();
      for (; p!= pEnd; p++) {
	speciesKey key = KeyConvert(p->second.dattrib[use_key]).getKey();
	iz = zval.find(key.first);
	if (iz != izEnd) speciesM[iz->second] += p->second.mass;
      }

      MPI_Allreduce(MPI_IN_PLACE, &speciesM[0], speciesM.size(), 
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      mol_weight = 0.0;
      j = 0;
      for (ZLtype::iterator i=ZList.begin(); i!=ZList.end(); i++) {
	mol_weight += speciesM[j++]/atomic_weights[*i];
      }
      mol_weight = 1.0/mol_weight;
    }
  }

  return mol_weight;
}
