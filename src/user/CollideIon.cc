#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cfloat>
#include <cmath>
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

// Warn if energy lost is smaller than COM energy available.  For
// debugging.  Set to false for production.
//
const bool frost_warning     = false;

// Very verbose selection debugging. Set to false for production.
//
const bool DEBUG_SL          = false;

// Verbose cross-section debugging. Set to false for production.
//
const bool DEBUG_CR          = false;

// Artifically suppress electron equipartition speed
//
const bool NO_DOF            = true;

// Artifically suppress electron equilibrium velocity
//
const bool NO_VEL            = false;

// Artifically prevent cooling by setting the energy removed from the
// COM frame to zero
//
const bool NO_COOL           = false;

// KE debugging; set to false for production
//
const bool KE_DEBUG          = true;

// Subtract KE from COM pair for testing only.  This is technically
// incorrect since the electrons are "trace" species and not part of
// the energy conservation.
//
const bool RECOMB_KE         = true;
const bool RECOMB_IP         = false;

// Cross-section debugging; set to false for production
//
const bool CROSS_DBG         = false;

// Excess trace map debugging; set to false for production
//
const bool EXCESS_DBG        = false;


// Minimum energy for Rutherford scattering of ions used to estimate
// the elastic scattering cross section
//
const double FloorEv = 0.05;

CollideIon::CollideIon(ExternalForce *force, Component *comp, 
		       double hD, double sD, 
		       const std::string& smap, int Nth) : 
  Collide(force, comp, hD, sD, Nth)
{
  // Read species file
  //
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
  dCross   .resize(nthrds);
  dInter   .resize(nthrds);
  sCross   .resize(nthrds);
  sInter   .resize(nthrds);
  meanF    .resize(nthrds);
  meanE    .resize(nthrds);
  meanR    .resize(nthrds);
  meanM    .resize(nthrds);
  neutF    .resize(nthrds);
  sCrsTot1 .resize(nthrds);
  sCrsTot2 .resize(nthrds);
  excessW  .resize(nthrds);
  CE1      .resize(nthrds);
  CE2      .resize(nthrds);
  kCE      .resize(nthrds);
  kEi      .resize(nthrds);
  kEe1     .resize(nthrds);
  kEe2     .resize(nthrds);
  Ein1     .resize(nthrds);
  Ein2     .resize(nthrds);
  Evel     .resize(nthrds);
  spTau    .resize(nthrds);
  spCrm    .resize(nthrds);
  spNsel   .resize(nthrds);
  spProb   .resize(nthrds);
  velER    .resize(nthrds);

  for (auto &v : velER) v.set_capacity(bufCap);

  //
  // Cross-section debugging [INIT]
  //
  if (CROSS_DBG) {
    nextTime_dbg = 0.0;		// Next target time
    nCnt_dbg     = 0;		// Number of cells accumulated so far

    ostringstream ostr;
    ostr << outdir << runtag << ".cross_section." << myid;
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
  labels[neut_neut  ] = "geometric ";
  labels[neut_elec  ] = "neutral el";
  labels[ion_elec   ] = "charged el";
  labels[free_free  ] = "free-free ";
  labels[colexcite  ] = "col excite";
  labels[ionize     ] = "ionization";
  labels[recomb     ] = "recombine ";

  labels[neut_neut_1] = "geometric  [1]";
  labels[neut_elec_1] = "neutral el [1]";
  labels[ion_elec_1 ] = "charged el [1]";
  labels[free_free_1] = "free-free  [1]";
  labels[colexcite_1] = "col excite [1]";
  labels[ionize_1   ] = "ionization [1]";
  labels[recomb_1   ] = "recombine  [1]";

  labels[neut_neut_2] = "geometric  [2]";
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
void CollideIon::initialize_cell(pHOT* const tree, pCell* const cell, 
				 double rvmax, int id)
{
				// Cache the calling tree
  curTree = tree;

  double KEtot, KEdspC;
  cell->KE(KEtot, KEdspC);	// KE in cell

  double massC = cell->Mass();	// Mass in cell

				// Add cell energy to diagnostic
				// handler
  collD->addCell(KEtot*massC, id);

				// Used for diagnostics only
  totalSoFar += massC * KEdspC;
  massSoFar  += massC;
  
  // Representative avg cell velocity in cgs
  //  
  double vavg = 0.5*rvmax*UserTreeDSMC::Vunit;

  vavg_dbg = 0.5*rvmax;

  // Representative avg cell energy in ergs
  //
  double Eerg = 0.5*vavg*vavg*amu;

  // In eV
  //
  double EeV = Eerg / eV;

  // Mean interparticle spacing in nm
  // 
  double ips = pow(cell->Volume()/cell->bods.size(), 0.333333) 
    * UserTreeDSMC::Lunit * 1.0e7;


  if (aType == Direct or aType == Weight) {

    unsigned eCnt  = 0;
    double   eVel1 = 0.0, eVel2 = 0.0;
    double   iVel1 = 0.0, iVel2 = 0.0;

    // Compute mean electron velocity
    //
    if (use_elec>=0) {
      
      // Sample cell
      //
      pCell *samp = cell->sample;

      if (samp) {
	for (auto c : samp->children) {
	  for (auto i : c.second->bods) {
	    Particle *p = c.second->Body(i);
	    KeyConvert k(p->iattrib[use_key]);

	    if (k.C()>1) {
	      for (int l=0; l<3; l++) {
		double ve  = p->dattrib[use_elec+l];
		eVel1 += ve;
		eVel2 += ve*ve;
		double vi  = p->vel[l];
		iVel1 += vi;
		iVel2 += vi*vi;
	      }
	      eCnt += k.C() - 1;
	    }
	  }
	}
      } else {
	for (auto i : cell->bods) {
	  Particle *p = cell->Body(i);
	  KeyConvert k(p->iattrib[use_key]);
	  
	  if (k.C()>1) {
	    for (int l=0; l<3; l++) {
	      double ve  = p->dattrib[use_elec+l];
	      eVel1 += ve;
	      eVel2 += ve*ve;
	      double vi  = p->vel[l];
	      iVel1 += vi;
	      iVel2 += vi*vi;
	    }
	    eCnt += k.C() - 1;
	  }
	}
      }

      if (eCnt>1) {
	eVel2 -= eVel1*eVel1/eCnt;
	iVel2 -= iVel1*iVel1/eCnt;
	Evel[id] = sqrt( fabs(eVel2 + iVel2)/(eCnt-1) );
      } else {
	Evel[id] = 0.0;
      }
    }

    // it1 and it2 are of type std::map<speciesKey, unsigned>; that
    // is, the number of particles of each speciesKey in the cell
    //
    for (auto it1 : cell->count) {

      speciesKey i1  = it1.first;
      double Radius1 = geometric(i1.first);
    
      // So, we are computing interactions for all possible
      // interaction pairs
      //
      for (auto it2 : cell->count) {
	
	speciesKey i2  = it2.first;
	double Radius2 = geometric(i2.first);

	double CrossG  = M_PI*(Radius1 + Radius2)*(Radius1 + Radius2);
	double Cross1  = 0.0;
	double Cross2  = 0.0;

	double m1      = atomic_weights[i1.first];
	double m2      = atomic_weights[i2.first];

	double ne1     = i1.second - 1;
	double ne2     = i2.second - 1;

	double dof1    = 1.0 + ne1;
	double dof2    = 1.0 + ne2;
	  
	if (NO_DOF) dof1 = dof2 = 1.0;

	double eVel1  = sqrt(atomic_weights[i1.first]/atomic_weights[0]/dof1);
	double eVel2  = sqrt(atomic_weights[i2.first]/atomic_weights[0]/dof2);

	if (use_elec) {
	  if (Evel[id]>0.0) eVel1 = eVel2 = Evel[id] / (0.5*rvmax);
	}

	if (NO_VEL) eVel1 = eVel2 = 1.0;

	if (i1.second>1 or i2.second>1) CrossG = 0.0;

	if (i2.second>1) {
	  if (i1.second==1)
	    Cross1 = elastic(i1.first, EeV*m1/dof2) * eVel2 * ne2;
	  else {
	    double b = 0.5*esu*esu*(i1.second - 1) /
	      std::max<double>(Eerg*m1/dof2, FloorEv*eV) * 1.0e7; // nm
	    b = std::min<double>(b, ips);
	    Cross1 = M_PI*b*b * eVel2 * ne2;
	  }
	}

	if (i1.second>1) {
	  if (i2.second==1)
	    Cross2 = elastic(i2.first, EeV*m2/dof1) * eVel1 * ne1;
	  else {
	    double b = 0.5*esu*esu*(i2.second - 1) /
	      std::max<double>(Eerg*m2/dof1, FloorEv*eV) * 1.0e7; // nm
	    b = std::min<double>(b, ips);
	    Cross2 = M_PI*b*b * eVel1 * ne1;
	  }
	}

	if (std::isnan(CrossG)) {
	  std::cout << "CrossG NaN" << std::endl;
	}

	if (std::isnan(Cross1)) {
	  std::cout << "Cross1 NaN" << std::endl;
	}

	if (std::isnan(Cross2)) {
	  std::cout << "Cross2 NaN" << std::endl;
	}

	csections[id][i1][i2] = (CrossG + Cross1 + Cross2) * crossfac * 1e-14 / 
	  (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
      }
    }
  }

  if (aType == Trace) {

    // In the trace computation, all superparticles are identical!
    // Hence, the use of the defaultKey.  This is slightly inefficient
    // but allows us to reuse the code for both the direct and trace
    // computations.
    //
    csections[id][defaultKey][defaultKey] = 0.0;

    // Compute mean weights in the cell
    //
    // In auto iterators below:
    //    s is of type std::map<speciesKey, int>
    //    b is of type std::vector<unsigned long>
    //
    // Per cell vvariables:
    //    meanF[id][sp] is the mean mass fraction for species sp
    //    meanE[id] is the mean number of electrons per particle
    //    meanR[id] is the mean effective cross-section radius
    //    neutF[id] is the neutral number fraction
    //    meanM[id] is the mean molecular weight
    //

				// Mean fraction in trace species
				// 
    for (auto s : SpList) meanF[id][s.first] = 0.0;
    neutF[id] = meanE[id] = meanR[id] = meanM[id] = 0.0;

				// Total mass of all particles and
				// relative fraction of trace species
				// and elctron number fraction in this
				// cell
    double massP = 0.0, numbP = 0.0;
    for (auto b : cell->bods) {
				// Particle mass accumulation
      Particle *p = tree->Body(b);
      massP      += p->mass;
				// Mass-weighted trace fraction
      for (auto s : SpList) {
	speciesKey k = s.first;
	double ee    = k.second - 1;
	double ww    = p->dattrib[s.second]/atomic_weights[k.first];
	
				// Mean mass fraction
	meanF[id][s.first] += p->mass * p->dattrib[s.second];
				// Mean electron number
	meanE[id]          += p->mass * ww * ee;
				// Mean number
	numbP              += p->mass * ww;

				// For neutrals only
	if (ee==0) {
				// Mean particle radius
	  meanR[id]        += p->mass * ww * geometric(k.first);
	  neutF[id]        += p->mass * ww;
	}
      }
    }
				// Normalize mass-weighted fraction
				// and electron fraction
				//
    if (massP>0.0) {
      for (auto s : SpList) meanF[id][s.first] /= massP;
      meanE[id] /= massP;
      neutF[id] /= massP;
    }
    if (neutF[id]>0.0) meanR[id] /= neutF[id];
    if (numbP    >0.0) meanM[id]  = massP/numbP;

				// Electron velocity factor for this
				// cell
    double eVel = sqrt(amu*meanM[id]/me);

				// Compute neutral and Coulombic cross
				// sections
    for (auto s : SpList) {

      speciesKey k = s.first;
				// Reduced mass for this interation
				// 
      double mu = 0.5*meanM[id];

				// Cumulative cross section
				// 
      double Cross = 0.0;

				// This species is a neutral
      if (k.second == 1) {
				// Neutral-neutral scattering 
				// cross section
	double Radius = geometric(k.first) + meanR[id];
	Cross += neutF[id] * M_PI*Radius*Radius;
	

				// Neutral-charged elastic cross section
	if (meanE[id] > 0.0)	// (recall Eerg and EeV are the mean 
				// interparticle KE)
	  Cross += elastic(k.first, EeV * mu) * eVel * meanE[id];

      } else {			// This species is an ion

				// Coulombic elastic scattering
	double b = 0.5*esu*esu*(k.second - 1) /
	  std::max<double>(Eerg*mu, FloorEv*eV) * 1.0e7; // nm
	b = std::min<double>(b, ips);
	Cross += M_PI*b*b * eVel * meanE[id];
      }

      double tCross = Cross * crossfac * 1e-14 / 
	(UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);

      csections[id][defaultKey][defaultKey] += tCross * 
	meanF[id][k] /atomic_weights[k.first];
    }
  }
}


sKey2Dmap& 
CollideIon::totalScatteringCrossSections(double crm, pCell* const c, int id)
{
  // it1 and it2 are of type std::map<speciesKey, unsigned>
  
  double vel  = crm * UserTreeDSMC::Vunit;
  double Eerg = 0.5*vel*vel*amu;
  double EeV  = Eerg / eV;

  // Mean interparticle spacing
  // 
  double ips = pow(c->Volume()/c->bods.size(), 0.333333) 
    * UserTreeDSMC::Lunit * 1.0e7;


  if (aType == Direct or aType == Weight) {

    for (auto it1 : c->count) {

      speciesKey i1 = it1.first;
      double geom1  = geometric(i1.first);
    
      for (auto it2 : c->count) {

	speciesKey i2 = it2.first;
	double geom2  = geometric(i2.first);
	  
	double m1     = atomic_weights[i1.first];
	double m2     = atomic_weights[i2.first];

	double ne1    = i1.second - 1;
	double ne2    = i2.second - 1;

	double dof1   = 1.0 + ne1;
	double dof2   = 1.0 + ne2;

	if (NO_DOF) dof1 = dof2 = 1.0;

	double eVel1  = sqrt(atomic_weights[i1.first]/atomic_weights[0]/dof1);
	double eVel2  = sqrt(atomic_weights[i2.first]/atomic_weights[0]/dof2);

	if (use_elec) eVel1 = eVel2 = Evel[id]/crm;

	if (NO_VEL)   eVel1 = eVel2 = 1.0;

	double Cross1 = 0.0;
	double Cross2 = 0.0;
	  
	// Both particles neutral?
	//
	if (i1.second==1 and i2.second==1) {
	  double Cross12 = M_PI*(geom1+geom2)*(geom1+geom2);
	  Cross1 = 0.5*Cross12;
	  Cross2 = 0.5*Cross12;
	}

	// Electrons in second particle?
	//
	if (ne2) {
	  if (i1.second==1)	// Neutral atom-electron scattering
	    Cross1 = elastic(i1.first, EeV*m1/dof2) * eVel2*ne2;
	  else {		// Rutherford scattering
	    double b = 0.5*esu*esu*(i1.second - 1) /
	      std::max<double>(Eerg*m1/dof2, FloorEv*eV) * 1.0e7; // nm
	    b = std::min<double>(b, ips);
	    Cross1 = M_PI*b*b * eVel2 * ne2;
	    }
	}

	// Electrons in first particle?
	//
	if (ne1) {
	    if (i2.second==1)	// Neutral atom-electron scattering
	      Cross2 = elastic(i2.first, EeV*m2/dof1) * eVel1*ne1;
	    else {		// Rutherford scattering
	      double b = 0.5*esu*esu*(i2.second - 1) /
		std::max<double>(Eerg*m2/dof1, FloorEv*eV) * 1.0e7; // nm
	      b = std::min<double>(b, ips);
	      Cross2 = M_PI*b*b * eVel1*ne1;
	    }
	}
	
	csections[id][i1][i2] = (Cross1 + Cross2) * crossfac * 1e-14 / 
	  (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
      }
    }

  }

  
  if (aType == Trace) {

    csections[id][defaultKey][defaultKey] = 0.0;

    // Compute the mean trace weight in the cell
    //
    // In auto iterators below:
    //    s is of type std::map<speciesKey, int>
    //    b is of type std::vector<unsigned long>
    //
				// Equipartition electron velocity
				//
    double eVel = sqrt(amu*meanM[id]/me);

				// Compute cross sections for all
				// interacting pairs
    for (auto s : SpList) {
      
      speciesKey k = s.first;
      double Cross = 0.0;

      if (k.second == 1) {
	// Default cross section for neutral cell
	//
	double Radius = geometric(k.first) + meanR[id];
	Cross += neutF[id] * M_PI*Radius*Radius;

	if (meanE[id]>0.0) {
      
	  // Use neutral or Rutherford scattering
	  //
	  Cross += elastic(k.first, EeV * meanM[id]) * eVel * meanE[id];
	}

      } else {			// This species is an ion
	
				// Coulumbic elastic scattering
	double b = 0.5*esu*esu*(k.second - 1) /
	  std::max<double>(Eerg*meanM[id], FloorEv*eV) * 1.0e7; // nm
	b = std::min<double>(b, ips);
	Cross += M_PI*b*b * eVel * meanE[id];
      }
      
      double tCross = Cross * crossfac * 1e-14 / 
	(UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
	
      csections[id][defaultKey][defaultKey] += tCross * 
	meanF[id][k]/atomic_weights[k.first];
    }
  }

  return csections[id];
}

double CollideIon::crossSectionDirect(pCell* const c, 
				      Particle* const p1, Particle* const p2, 
				      double cr, int id)
{
  // Species keys
  //
  KeyConvert k1(p1->iattrib[use_key]);
  KeyConvert k2(p2->iattrib[use_key]);

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

  // Available COM energy
  //
  kEi[id] = 0.5 * mu * vel*vel;

  // Electron velocity equipartition factors
  //
  double eVel = sqrt(mu/me);

  if (NO_VEL) eVel = 1.0;

  // Mean interparticle spacing
  // 
  double ips = pow(c->Volume()/c->bods.size(), 0.333333) 
    * UserTreeDSMC::Lunit * 1.0e7;

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
  
  // Save the per-interaction cross sections
  //
  dCross[id].clear();

  // Index the interactions
  //
  dInter[id].clear();

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
    dCross[id].push_back(cross12*crossfac);
    dInter[id].push_back(neut_neut_1);

    cross21 = geometric(Z2);
    dCross[id].push_back(cross21*crossfac);
    dInter[id].push_back(neut_neut_2);
  }

				//-------------------------------
				// Electrons in second particle
				//-------------------------------
  if (ne2 > 0) {
    if (C1==1) {		// Neutral atom-electron scattering
      cross12 = elastic(Z1, kEe2[id]) * eVel*ne2 * crossfac;
      dCross[id].push_back(cross12);
      dInter[id].push_back(neut_elec_1);
    }  else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C1-1) /
	std::max<double>(kEe2[id]*eV, FloorEv*eV) * 1.0e7; // nm
      b = std::min<double>(b, ips);
      cross12 = M_PI*b*b * eVel*ne2 * crossfac;
      dCross[id].push_back(cross12);
      dInter[id].push_back(ion_elec_1);
    }
  }
    
				//-------------------------------
				// Electrons in first particle
				//-------------------------------
  if (ne1 > 0) {
    if (C2==1) {		// Neutral atom-electron scattering
      cross21 = elastic(Z2, kEe1[id]) * eVel*ne1 * crossfac;
      dCross[id].push_back(cross21);
      dInter[id].push_back(neut_elec_2);
    } else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C2-1) /
	std::max<double>(kEe1[id]*eV, FloorEv*eV) * 1.0e7; // nm
      b = std::min<double>(b, ips);
      cross21 = M_PI*b*b * eVel*ne1 * crossfac;
      dCross[id].push_back(cross21);
      dInter[id].push_back(ion_elec_2);
    }
  }

  //--------------------------------------------------
  // Ion keys
  //--------------------------------------------------

  lQ Q1(Z1, C1), Q2(Z2, C2);

  //===================================================================
  //  ___      _                                      _   _    _     
  // | _ \_  _| |_   _ _  _____ __ __  _ __  __ _ _ _| |_(_)__| |___ 
  // |  _/ || |  _| | ' \/ -_) V  V / | '_ \/ _` | '_|  _| / _| / -_)
  // |_|  \_,_|\__| |_||_\___|\_/\_/  | .__/\__,_|_|  \__|_\__|_\___|
  //                                  |_|                            
  //  _     _                   _   _               _                
  // (_)_ _| |_ ___ _ _ __ _ __| |_(_)___ _ _  ___ | |_  ___ _ _ ___ 
  // | | ' \  _/ -_) '_/ _` / _|  _| / _ \ ' \(_-< | ' \/ -_) '_/ -_)
  // |_|_||_\__\___|_| \__,_\__|\__|_\___/_||_/__/ |_||_\___|_| \___|
  //                                                                 
  //===================================================================
  

  //--------------------------------------------------
  // Particle 1 interacts with Particle 2
  //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
  if (C1 > 1 and ne2 > 0) {	// Ion and Ion only

    double ff1 = ch.IonList[Q1]->freeFreeCross(kEe2[id], id);
    double crs = eVel*ne2 * ff1;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(free_free_1);
      sum12 += crs;
    }
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  if (ne2 > 0 and C1 <= Z1) {	// Particle 1 must be bound

    CE1[id] = ch.IonList[Q1]->collExciteCross(kEe2[id], id);

    double crs = eVel*ne2 * CE1[id].back().first;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(colexcite_1);
      sum12 += crs;
    }
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  if (ne2 > 0 and C1 <= Z1) {	// Particle 1 must be bound

    double DI1 = ch.IonList[Q1]->directIonCross(kEe2[id], id);
    double crs = eVel*ne2 * DI1;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(ionize_1);
      sum12 += crs;
    }
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  if (C1 > 1 and ne2 > 0) {	// Particle 1 must be an ion

    std::vector<double> RE1 = ch.IonList[Q1]->radRecombCross(kEe2[id], id);
    double crs = eVel*ne2 * RE1.back();

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(recomb_1);
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
    double ff2 = ch.IonList[Q2]->freeFreeCross(kEe1[id], id);
    double crs = eVel*ne1 * ff2;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(free_free_2);
      sum21 += crs;
    }
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  if (ne1 > 0 and C2 <= Z2) {

    CE2[id] = ch.IonList[Q2]->collExciteCross(kEe1[id], id);
    double crs = eVel*ne1 * CE2[id].back().first;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(colexcite_2);
      sum21 += crs;
    }
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  if (ne1 > 0 and C2 <= Z2) {
    double DI2 = ch.IonList[Q2]->directIonCross(kEe1[id], id);
    double crs = ne1 * DI2;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(ionize_2);
      sum21 += crs;
    }
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  if (C2 > 1 and ne1 > 0) {
    std::vector<double> RE2 = ch.IonList[Q2]->radRecombCross(kEe1[id], id);
    double crs = eVel*ne1*RE2.back();

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(recomb_2);
      sum21 += crs;
    } 
  }

				//-------------------------------
				// *** Convert to system units
				//-------------------------------
  return (cross12 + cross21 + sum12 + sum21) * 1e-14 / 
    (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
}

double CollideIon::crossSectionWeight(pCell* const c, 
				      Particle* const _p1, Particle* const _p2, 
				      double cr, int id)
{
  Particle* p1 = _p1;		// Pointer copies
  Particle* p2 = _p2;

  // Mean interparticle spacing
  // 
  double ips = pow(c->Volume()/c->bods.size(), 0.333333) 
    * UserTreeDSMC::Lunit * 1.0e7;

  // Species keys
  //
  KeyConvert k1(p1->iattrib[use_key]);
  KeyConvert k2(p2->iattrib[use_key]);

  unsigned short Z1 = k1.getKey().first, C1 = k1.getKey().second;
  unsigned short Z2 = k2.getKey().first, C2 = k2.getKey().second;
  
  if (ZWList[Z1] < ZWList[Z2]) {
    Particle *pT = p1;
    p1 = p2;
    p2 = pT;

    k1 = KeyConvert(p1->iattrib[use_key]);
    k2 = KeyConvert(p2->iattrib[use_key]);

    Z1 = k1.getKey().first;
    C1 = k1.getKey().second;
    Z2 = k2.getKey().first;
    C2 = k2.getKey().second;
  }
      
  // Number of atoms in each super particle
  //
  double N1   = p1->mass*UserTreeDSMC::Munit/amu / atomic_weights[Z1];
  double N2   = p2->mass*UserTreeDSMC::Munit/amu / atomic_weights[Z2];

  // Number of associated electrons for each particle
  //
  double ne1  = C1 - 1;
  double ne2  = C2 - 1;
  
  double dof1 = 1.0 + ne1;
  double dof2 = 1.0 + ne2;

  if (NO_DOF) dof1 = dof2 = 1.0;

  // Energy available in the center of mass of the atomic collision
  //
  double vel = cr * UserTreeDSMC::Vunit;
  double m1  = atomic_weights[Z1]*amu;
  double m2  = atomic_weights[Z2]*amu;
  double me  = atomic_weights[ 0]*amu;
  double mu0 = m1 * m2 / (m1 + m2);
  double mu1 = m1;
  double mu2 = m2;

  // Electron velocity equipartition factors
  //
  double eVel1 = sqrt(m1/me/dof1);
  double eVel2 = sqrt(m2/me/dof2);

  if (NO_VEL) {
    eVel1 = eVel2 = 1.0;
  } else if (use_elec) {
    eVel1 = eVel2 = 0.0;
    for (unsigned i=0; i<3; i++) {
      double rvel1 = p1->dattrib[use_elec+i] - p2->vel[i];
      double rvel2 = p2->dattrib[use_elec+i] - p1->vel[i];
      eVel1 += rvel1*rvel1;
      eVel2 += rvel2*rvel2;
    }
    eVel1 = sqrt(eVel1) * UserTreeDSMC::Vunit;
    eVel2 = sqrt(eVel2) * UserTreeDSMC::Vunit;
  }

  // Available COM energy
  //
  kEi[id] = 0.5 * mu0 * vel*vel;

  if (use_elec) {
    kEe1[id] = 0.5 * me * eVel2*eVel2;
    kEe2[id] = 0.5 * me * eVel1*eVel1;
  } else {
    kEe1[id] = 0.5 * mu1 * vel*vel/dof2;
    kEe2[id] = 0.5 * mu2 * vel*vel/dof1;
  }

  // These are now ratios
  //
  eVel1 /= vel;
  eVel2 /= vel;

  // Internal energy per particle
  //
  Ein1[id] = Ein2[id] = 0.0;

  if (use_Eint>=0) {
    Ein1[id] = p1->dattrib[use_Eint] * UserTreeDSMC::Eunit / N1;
    Ein2[id] = p2->dattrib[use_Eint] * UserTreeDSMC::Eunit / N2;

    // Compute the total available energy and divide among degrees of freedom
    // Convert ergs to eV
    //
    kEe1[id] = (kEe1[id] + Ein1[id]) / eV;
    kEe2[id] = (kEe1[id] + Ein2[id]) / eV;
  } else {
    kEe1[id] /= eV;
    kEe2[id] /= eV;
  }
  
  kEi[id] /= eV;

  // Save the per-interaction cross sections
  //
  dCross[id].clear();

  // Index the interactions
  //
  dInter[id].clear();

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
    dCross[id].push_back(cross12*crossfac);
    dInter[id].push_back(neut_neut_1);

    cross21 = geometric(Z2);
    dCross[id].push_back(cross21*crossfac);
    dInter[id].push_back(neut_neut_2);
  }

				//-------------------------------
				// Electrons in second particle
				//-------------------------------
  if (ne2 > 0) {
    if (C1==1) {		// Neutral atom-electron scattering
      cross12 = elastic(Z1, kEe1[id]) * eVel2 * ne2 * crossfac;
      dCross[id].push_back(cross12);
      dInter[id].push_back(neut_elec_1);
    }  else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C1-1) /
	std::max<double>(kEe1[id]*eV, FloorEv*eV) * 1.0e7; // nm
      b = std::min<double>(b, ips);
      cross12 = M_PI*b*b * eVel2 * ne2 * crossfac;
      dCross[id].push_back(cross12);
      dInter[id].push_back(ion_elec_1);
    }
  }
    
				//-------------------------------
				// Electrons in first particle
				//-------------------------------
  if (ne1 > 0) {
    if (C2==1) {		// Neutral atom-electron scattering
      cross21 = elastic(Z2, kEe2[id]) * eVel1 * ne1 * crossfac;
      dCross[id].push_back(cross21);
      dInter[id].push_back(neut_elec_2);
    } else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C2-1) /
	std::max<double>(kEe2[id]*eV, FloorEv*eV) * 1.0e7; // nm
      b = std::min<double>(b, ips);
      cross21 = M_PI*b*b * eVel1 * ne1 * crossfac;
      dCross[id].push_back(cross21);
      dInter[id].push_back(ion_elec_2);
    }
  }

  //--------------------------------------------------
  // Ion keys
  //--------------------------------------------------

  lQ Q1(Z1, C1), Q2(Z2, C2);

  //===================================================================
  //  ___      _                                      _   _    _     
  // | _ \_  _| |_   _ _  _____ __ __  _ __  __ _ _ _| |_(_)__| |___ 
  // |  _/ || |  _| | ' \/ -_) V  V / | '_ \/ _` | '_|  _| / _| / -_)
  // |_|  \_,_|\__| |_||_\___|\_/\_/  | .__/\__,_|_|  \__|_\__|_\___|
  //                                  |_|                            
  //  _     _                   _   _               _                
  // (_)_ _| |_ ___ _ _ __ _ __| |_(_)___ _ _  ___ | |_  ___ _ _ ___ 
  // | | ' \  _/ -_) '_/ _` / _|  _| / _ \ ' \(_-< | ' \/ -_) '_/ -_)
  // |_|_||_\__\___|_| \__,_\__|\__|_\___/_||_/__/ |_||_\___|_| \___|
  //                                                                 
  //===================================================================
  

  //--------------------------------------------------
  // Particle 1 interacts with Particle 2
  //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
  if (C1 > 1 and ne2 > 0) {	// Ion and Ion only

    double ff1 = ch.IonList[Q1]->freeFreeCross(kEe1[id], id);
    double crs = eVel2 * ne2 * ff1;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(free_free_1);
      sum12 += crs;
    }
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  if (ne2 > 0 and C1 <= Z1) {	// Particle 1 must be bound

    CE1[id] = ch.IonList[Q1]->collExciteCross(kEe1[id], id);

    double crs = eVel2 * ne2 * CE1[id].back().first;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(colexcite_1);
      sum12 += crs;
    }
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  if (ne2 > 0 and C1 <= Z1) {	// Particle 1 must be bound

    double DI1 = ch.IonList[Q1]->directIonCross(kEe1[id], id);
    double crs = eVel2 * ne2 * DI1;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(ionize_1);
      sum12 += crs;
    }
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  if (C1 > 1 and ne2 > 0) {	// Particle 1 must be an ion

    std::vector<double> RE1 = ch.IonList[Q1]->radRecombCross(kEe1[id], id);
    double crs = eVel2 * ne2 * RE1.back();

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(recomb_1);
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
    double ff2 = ch.IonList[Q2]->freeFreeCross(kEe2[id], id);
    double crs = eVel1 * ne1 * ff2;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(free_free_2);
      sum21 += crs;
    }
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  if (ne1 > 0 and C2 <= Z2) {

    CE2[id] = ch.IonList[Q2]->collExciteCross(kEe2[id], id);
    double crs = eVel1 * ne1 * CE2[id].back().first;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(colexcite_2);
      sum21 += crs;
    }
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  if (ne1 > 0 and C2 <= Z2) {
    double DI2 = ch.IonList[Q2]->directIonCross(kEe2[id], id);
    double crs = eVel1 * ne1 * DI2;

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(ionize_2);
      sum21 += crs;
    }
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  if (C2 > 1 and ne1 > 0) {
    std::vector<double> RE2 = ch.IonList[Q2]->radRecombCross(kEe2[id], id);
    double crs = eVel1 * ne1 * RE2.back();

    if (crs>0.0) {
      dCross[id].push_back(crs);
      dInter[id].push_back(recomb_2);
      sum21 += crs;
    } 
  }

				//-------------------------------
				// *** Convert to system units
				//-------------------------------
  return (cross12 + cross21 + sum12 + sum21) * 1e-14 / 
    (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
}


double CollideIon::crossSectionTrace(pCell* const c, 
				     Particle* const p1, Particle* const p2, 
				     double cr, int id)
{
  double totalCross = 0.0;

  //
  // Clear excess species change map
  //
  excessW[id].clear();

  //
  // Clear total cross section map
  //
  sCrsTot1[id].clear();
  sCrsTot2[id].clear();

  //
  // Compute "matrix" of interactions
  //
  
  // s1 and s2 are of type std::map<speciesKey, int>

  // Mean interparticle spacing
  // 
  double ips = pow(c->Volume()/c->bods.size(), 0.333333) 
    * UserTreeDSMC::Lunit * 1.0e7;

  // Translational COM energy
  //
  double vel = cr * UserTreeDSMC::Vunit;

  kEi[id] = 0.5 * meanM[id] * vel*vel * amu;

  // Convert the total available energy from ergs to eV
  //
  kEe1[id] = kEi[id] / eV;
  kEe2[id] = kEi[id] / eV;

  double kEe = kEi[id] / eV;

  // Electron velocity equipartition factors
  //
  double eVel = sqrt(meanM[id]*amu/me);

  for (auto s : SpList) {

    // Particle species: key and number weight
    //
    speciesKey k = s.first;
    
    // Atomic numbers
    //
    unsigned short Z = k.first, C = k.second;
  
    // Save the per-interaction cross sections
    std::vector<double> tCrossMap;
      
    // Index the interactions
    std::vector<int> tInterMap;

    //--------------------------------------------------
    // Total scattering cross section
    //--------------------------------------------------
    
    double crossS = 0.0;

    //--------------------------------------------------
    // Total inelastic cross section
    //--------------------------------------------------
    
    double crossI = 0.0;

    //--------------------------------------------------
    // Number weighting for this particle
    //--------------------------------------------------
    
    double w1 = p1->dattrib[SpList[k]];
    double w2 = p2->dattrib[SpList[k]];
    double ww = (w1 + w2)/atomic_weights[k.first];

    //-------------------------------
    // Elastic: particle is neutral
    //-------------------------------

    if (C==1) {

      //
      // Neutral-neutral geometric cross section
      //

				// Geometric cross sections based on
				// atomic radius
      double Radius = geometric(Z) + meanR[id];

      crossS += neutF[id] * M_PI*Radius*Radius * crossfac;
      tCrossMap.push_back(crossS);
      tInterMap.push_back(neut_neut);
	
      //
      // Neutral atom-electron scattering
      //
      
      crossS += elastic(Z, kEe) * eVel * meanE[id] * crossfac;
      
      tCrossMap.push_back(crossS);
      tInterMap.push_back(neut_elec);
      
    } else {
				// 
				// Rutherford scattering
				//
      double b = 0.5*esu*esu*(C-1) /
	std::max<double>(kEe*eV, FloorEv*eV) * 1.0e7; // nm
      b = std::min<double>(b, ips);
      crossS += M_PI*b*b * eVel * meanE[id] * crossfac;

      tCrossMap.push_back(crossS);
      tInterMap.push_back(ion_elec);
    }

    //--------------------------------------------------
    // Ion key
    //--------------------------------------------------
    lQ Q(Z, C);


    //===================================================================
    //  ___      _                                      _   _    _     
    // | _ \_  _| |_   _ _  _____ __ __  _ __  __ _ _ _| |_(_)__| |___ 
    // |  _/ || |  _| | ' \/ -_) V  V / | '_ \/ _` | '_|  _| / _| / -_)
    // |_|  \_,_|\__| |_||_\___|\_/\_/  | .__/\__,_|_|  \__|_\__|_\___|
    //                                  |_|                            
    //  _     _                   _   _               _                
    // (_)_ _| |_ ___ _ _ __ _ __| |_(_)___ _ _  ___ | |_  ___ _ _ ___ 
    // | | ' \  _/ -_) '_/ _` / _|  _| / _ \ ' \(_-< | ' \/ -_) '_/ -_)
    // |_|_||_\__\___|_| \__,_\__|\__|_\___/_||_/__/ |_||_\___|_| \___|
    //                                                                 
    //===================================================================

				//-------------------------------
				// *** Free-free
				//-------------------------------
    if (C > 1 and meanE[id] > 0) {
      double ff = ch.IonList[Q]->freeFreeCross(kEe, id);
      double crs = meanE[id] * eVel * ff;

      if (crs>0.0) {
	tCrossMap.push_back(crs);
	tInterMap.push_back(free_free);
	crossI += crs;
      }
    }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
				// Particle must be bound
    if (meanE[id] > 0 and C <= Z) {

      CEvector V = ch.IonList[Q]->collExciteCross(kEe, id);
      double crs = meanE[id] * eVel * V.back().first;

      if (crs>0.0) {
	tCrossMap.push_back(crs);
	tInterMap.push_back(colexcite);
	kCE[id][k] = V;
	crossI += crs;
      }
    }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
				// Particle must be bound
    if (meanE[id] > 0 and C <= Z) {

      double DI = ch.IonList[Q]->directIonCross(kEe, id);
      double crs = meanE[id] * eVel * DI;

      if (crs>0.0) {
	tCrossMap.push_back(crs);
	tInterMap.push_back(ionize);
	crossI += crs;
      }
    }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
				// Particle must be an ion
    if (C > 1 and meanE[id] > 0) {

      std::vector<double> RE = ch.IonList[Q]->radRecombCross(kEe, id);
      double crs = meanE[id] * eVel * RE.back();

      if (crs>0.0) {
	tCrossMap.push_back(crs);
	tInterMap.push_back(recomb);
	crossI += crs;
      }
    }

    sCross[id][k] = tCrossMap;
    sInter[id][k] = tInterMap;

    totalCross += (crossS + crossI)*ww;
  }
  
  return totalCross;
}


int CollideIon::inelasticDirect(pCell* const c, 
				Particle* const p1, Particle* const p2, 
				double *cr, int id)
{
  int ret = 0;			// No error (flag)
  int interFlag = -1;		// Invalid value by default

  // Species keys
  //
  KeyConvert k1(p1->iattrib[use_key]);
  KeyConvert k2(p2->iattrib[use_key]);

  collTDPtr ctd1 = (*collD)[k1.getKey()];
  collTDPtr ctd2 = (*collD)[k2.getKey()];

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
  double dE    = kE*TolV*TolV;
  double delE  = 0.0;

  // Now that the interactions have been calculated, create the
  // normalized cross section list to pick the interaction
  //
  std::vector<double> TotalCross;
  double tCross = 0.0;
  int si = 0;
  for (size_t i = 0; i < dCross[id].size(); i++) {
    if (std::isnan(dCross[id][i])) {
      std::cout << "dCross[" << id << "][" << i << "] is Nan"
		<< std::setw(14) << dInter[id][i]
		<< std::setw(18) << labels[dInter[id][i]]
		<< std::endl;
    } else if (std::isinf(dCross[id][i])) {
      std::cout << "dCross[" << id << "][" << i << "] is "
		<< dCross[id][i]
		<< std::setw(14) << dInter[id][i]
		<< std::setw(18) << labels[dInter[id][i]]
		<< std::endl;
    } else {
      tCross += dCross[id][i];
      TotalCross.push_back(tCross);
      si++;
    }
  }

  int partflag = 0;		// Will be 1 or 2, dependending on
				// which ion or neutral is selected
				// for interaction.  Will be 0 if no
				// interaction is selected.
  
  if (tCross != 0) {
    // Cumulative cross-section distribution for interaction selection
    //
    std::vector<double> CDF;
    for (size_t i = 0; i < TotalCross.size(); i++) {
      if (std::isnan(TotalCross[i])) {
	std::cout << "TotalCross[i][" << id << "][" << i << "] is Nan"
		  << std::endl;
      } else {
	CDF.push_back(TotalCross[i]/tCross);
      }
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
		<< ", ran=" << ran 
		<< ", siz=" << CDF.size()
		<< ", beg=" << CDF.front()
		<< ", end=" << CDF.back()
		<< ", tot=" << tCross
		<< std::endl;
      index = 0;
    }

    // Finally, set the interaction type based on the selected index
    //
    interFlag = dInter[id][index];

    //-------------------------
    // VERBOSE DEBUG TEST
    //-------------------------
    // Set to false for production
    //  |
    //  v
    if (DEBUG_CR) {
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
	for (size_t i = 0; i < dCross[id].size(); i++) {
	  std::cout << std::setw( 8) << i
		    << std::setw( 8) << dInter[id][i]
		    << std::setw(14) << dCross[id][i]
		    << std::setw(14) << CDF[i]
		    << std::setw(18) << labels[dInter[id][i]]
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
      delE          = IS.selectFFInteract(ch.IonList[Q1], id);
      partflag      = 1;
      std::get<0>(ctd1->ff[id])++; 
      std::get<1>(ctd1->ff[id]) += NN;
      std::get<2>(ctd1->ff[id]) += delE * NN;
    }

    if (interFlag == colexcite_1) {
      delE = IS.selectCEInteract(ch.IonList[Q1], CE1[id]);
      partflag      = 1;
      std::get<0>(ctd1->CE[id])++;
      std::get<1>(ctd1->CE[id]) += NN;
      std::get<2>(ctd1->CE[id]) += delE * NN;
    }

    if (interFlag == ionize_1) {
      delE          = IS.DIInterLoss(ch.IonList[Q1]);
      p1->iattrib[use_key] = k1.updateC(++C1);
      partflag      = 1;
      std::get<0>(ctd1->CI[id])++; 
      std::get<1>(ctd1->CI[id]) += NN;
      std::get<2>(ctd1->CI[id]) += delE * NN;
    }

    if (interFlag == recomb_1) {
      if (RECOMB_KE) {
	if (RECOMB_IP) {
	  lQ rQ(Z1, C1-1);
	  delE = ch.IonList[rQ]->ip + kEe1[id];
	} else 
	  delE = kEe1[id];
      }

      p1->iattrib[use_key] = k1.updateC(--C1);
      partflag      = 1;
      std::get<0>(ctd1->RR[id])++; 
      std::get<1>(ctd1->RR[id]) += NN;
      std::get<2>(ctd1->RR[id]) += delE * NN;
    }
    
    //-------------------------
    // Particle 2 interactions
    //-------------------------

    if (interFlag == free_free_2) {
      delE          = IS.selectFFInteract(ch.IonList[Q2], id);
      partflag      = 2;
      std::get<0>(ctd2->ff[id])++;
      std::get<1>(ctd2->ff[id]) += NN;
      std::get<2>(ctd2->ff[id]) += delE * NN;
    }

    if (interFlag == colexcite_2) {
      delE         = IS.selectCEInteract(ch.IonList[Q2], CE2[id]);
      partflag     = 2;
      std::get<0>(ctd2->CE[id])++; 
      std::get<1>(ctd2->CE[id]) += NN;
      std::get<2>(ctd2->CE[id]) += delE * NN;
    }

    if (interFlag == ionize_2) {
      delE = IS.DIInterLoss(ch.IonList[Q2]);
      p2->iattrib[use_key] = k2.updateC(++C2);
      std::get<0>(ctd2->CI[id])++; 
      std::get<1>(ctd2->CI[id]) += NN;
      std::get<2>(ctd2->CI[id]) += delE * NN;
      partflag     = 2;
    }

    if (interFlag == recomb_2) {
      if (RECOMB_KE) {
	if (RECOMB_IP) {
	  lQ rQ(Z2, C2-1);
	  delE = ch.IonList[rQ]->ip + kEe2[id];
	} else
	  delE = kEe2[id];
      }

      p2->iattrib[use_key] = k2.updateC(--C2);
      partflag     = 2;
      std::get<0>(ctd2->RR[id])++; 
      std::get<1>(ctd2->RR[id]) += NN;
      std::get<2>(ctd2->RR[id]) += delE * NN;
    }

    // Convert to super particle
    //
    if (partflag) delE *= NN;
    
    // Convert back to cgs
    //
    delE = delE * eV;
  }
  
  // For elastic interactions, delE == 0
  //
  assert(delE >= 0.0);

  // Artifically prevent cooling by setting the energy removed from
  // the COM frame to zero
  //
  if (NO_COOL) delE = 0.0;

  // Elastic event
  //
  if (delE<=0.0) return ret;

  // Convert energy loss to system units
  //
  delE = delE/UserTreeDSMC::Eunit;
  
  // Assign interaction energy variables
  //
  double remE=0.0, totE=0.0, kEe=0.0;

  // Electrons from Particle 2 have interacted with atom/ion in Particle 1
  //
  if (partflag==1) {
    totE = kE;			// KE + internal
    if (use_Eint>=0) totE += p2->dattrib[use_Eint];
    remE = totE - dE;		// Energy floor
    kEe  = kEe2[id];		// Electron energy

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
    double dtE = remE/delE * spTau[id];
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
    double kEm = totE;
    
				// Get new relative velocity
    (*cr)          = sqrt( 2.0*kEm/Mu );

    ret            = 0;		// No error

    if (partflag==1) {
      std::get<0>(ctd1->dv[id])++; 
      std::get<1>(ctd1->dv[id]) += N1;
      std::get<2>(ctd1->dv[id]) += 
	0.5*Mu*(vi - (*cr))*(vi - (*cr))/N1 * UserTreeDSMC::Eunit / eV;
    }
    
    if (partflag==2) {
      std::get<0>(ctd2->dv[id])++; 
      std::get<1>(ctd2->dv[id]) += N2;
      std::get<2>(ctd2->dv[id]) += 
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
	std::get<0>(ctd1->dv[id])++; 
	std::get<1>(ctd1->dv[id]) += N1;
	std::get<2>(ctd1->dv[id]) += 
	  0.5*Mu*(vi - (*cr))*(vi - (*cr))/N1 * UserTreeDSMC::Eunit / eV;
      }
      
      if (partflag==2) {
	std::get<0>(ctd2->dv[id])++; 
	std::get<1>(ctd2->dv[id]) += N2; 
	std::get<2>(ctd2->dv[id]) +=
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
	std::get<0>(ctd1->dv[id])++; 
	std::get<1>(ctd1->dv[id]) += N1;
	std::get<2>(ctd1->dv[id]) +=
	  0.5*Mu*(vi - (*cr))*(vi - (*cr))/N1 * UserTreeDSMC::Eunit / eV;
      }
      
      if (partflag==2) {
	std::get<0>(ctd2->dv[id])++; 
	std::get<1>(ctd2->dv[id]) += N2;
	std::get<2>(ctd2->dv[id]) +=
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

int CollideIon::inelasticWeight(pCell* const c, 
				Particle* const _p1, Particle* const _p2,
				double *cr, int id)
{
  int ret = 0;			// No error (flag)
  int interFlag = -1;		// Invalid value by default

  Particle* p1 = _p1;
  Particle* p2 = _p2;

  // Species keys
  //
  KeyConvert k1(p1->iattrib[use_key]);
  KeyConvert k2(p2->iattrib[use_key]);

  collTDPtr ctd1 = (*collD)[k1.getKey()];
  collTDPtr ctd2 = (*collD)[k2.getKey()];

  unsigned short Z1 = k1.getKey().first, C1 = k1.getKey().second;
  unsigned short Z2 = k2.getKey().first, C2 = k2.getKey().second;

  if (p1->mass/atomic_weights[Z1] < p2->mass/atomic_weights[Z2]) {

    // Swap the particle pointers
    //
    Particle *pT = p1;
    p1 = p2;
    p2 = pT;

    // Swap the collision diag pointers
    //
    collTDPtr ctdT = ctd1;
    ctd1 = ctd2;
    ctd2 = ctdT;

    // Reassign the keys and species indices
    //
    k1 = KeyConvert(p1->iattrib[use_key]);
    k2 = KeyConvert(p2->iattrib[use_key]);

    Z1 = k1.getKey().first;
    C1 = k1.getKey().second;

    Z2 = k2.getKey().first;
    C2 = k2.getKey().second;
  }
      
  // Find the trace ratio
  //
  double Wa = p1->mass / atomic_weights[Z1];
  double Wb = p2->mass / atomic_weights[Z2];
  double  q = Wb / Wa;

  // Number interacting atoms
  //
  double NN = Wb * UserTreeDSMC::Munit / amu;

  // For tracking energy conservation (system units)
  //
  double delE  = 0.0;

  // Now that the interactions have been calculated, create the
  // normalized cross section list to pick the interaction
  //
  std::vector<double> TotalCross;
  double tCross = 0.0;
  int si = 0;
  for (size_t i = 0; i < dCross[id].size(); i++) {
    if (std::isnan(dCross[id][i])) {
      std::cout << "dCross[" << id << "][" << i << "] is Nan"
		<< std::setw(14) << dInter[id][i]
		<< std::setw(18) << labels[dInter[id][i]]
		<< std::endl;
    } else if (std::isinf(dCross[id][i])) {
      std::cout << "dCross[" << id << "][" << i << "] is "
		<< dCross[id][i]
		<< std::setw(14) << dInter[id][i]
		<< std::setw(18) << labels[dInter[id][i]]
		<< std::endl;
    } else {
      tCross += dCross[id][i];
      TotalCross.push_back(tCross);
      si++;
    }
  }

  int partflag = 0;		// Will be 1 or 2, dependending on
				// which ion or neutral is selected
				// for interaction.  Will be 0 if no
				// interaction is selected.

  if (tCross != 0) {
    // Cumulative cross-section distribution for interaction selection
    //
    std::vector<double> CDF;
    for (size_t i = 0; i < TotalCross.size(); i++) {
      if (std::isnan(TotalCross[i])) {
	std::cout << "TotalCross[i][" << id << "][" << i << "] is Nan"
		  << std::endl;
      } else {
	CDF.push_back(TotalCross[i]/tCross);
      }
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
		<< ", ran=" << ran 
		<< ", siz=" << CDF.size()
		<< ", beg=" << CDF.front()
		<< ", end=" << CDF.back()
		<< ", tot=" << tCross
		<< std::endl;
      index = 0;
    }

    // Finally, set the interaction type based on the selected index
    //
    interFlag = dInter[id][index];

    //-------------------------
    // VERBOSE DEBUG TEST
    //-------------------------
    // Set to false for production
    //  |
    //  v
    if (DEBUG_CR) {
      speciesKey i1 = k1.getKey();
      speciesKey i2 = k2.getKey();
      double cfac = 1e-14 / (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
      //
      // Output on collisions for now . . . 
      //
      std::cout << std::setw( 8) << "index"
		<< std::setw( 8) << "flag"
		<< std::setw(14) << "cross"
		<< std::setw(14) << "cumul"
		<< std::setw(14) << "tCross"
		<< std::setw(18) << "type label"
		<< std::endl
		<< std::setw( 8) << "-----"
		<< std::setw( 8) << "-----"
		<< std::setw(14) << "---------"
		<< std::setw(14) << "---------"
		<< std::setw(14) << "---------"
		<< std::setw(18) << "---------------"
		<< std::endl;
      for (size_t i = 0; i < dCross[id].size(); i++) {
	std::cout << std::setw( 8) << i
		  << std::setw( 8) << dInter[id][i]
		  << std::setw(14) << dCross[id][i]
		  << std::setw(14) << CDF[i]
		  << std::setw(14) << TotalCross[i]/csections[id][i1][i2] * cfac
		  << std::setw(18) << labels[dInter[id][i]]
		  << std::endl;
      }
      std::cout << std::endl;

    }

    //--------------------------------------------------
    // Ion keys
    //--------------------------------------------------

    lQ Q1(Z1, C1), Q2(Z2, C2);

    //-------------------------
    // Particle 1 interactions
    //-------------------------

    if (interFlag == free_free_1) {
      delE          = IS.selectFFInteract(ch.IonList[Q1], id);
      partflag      = 1;
      std::get<0>(ctd1->ff[id])++; 
      std::get<1>(ctd1->ff[id]) += Wb;
      std::get<2>(ctd1->ff[id]) += delE * NN;
    }

    if (interFlag == colexcite_1) {
      delE = IS.selectCEInteract(ch.IonList[Q1], CE1[id]);
      partflag      = 1;
      std::get<0>(ctd1->CE[id])++;
      std::get<1>(ctd1->CE[id]) += Wb;
      std::get<2>(ctd1->CE[id]) += delE * NN;
    }

    if (interFlag == ionize_1) {
      delE          = IS.DIInterLoss(ch.IonList[Q1]);
      p1->iattrib[use_key] = k1.updateC(++C1);
      partflag      = 1;
      std::get<0>(ctd1->CI[id])++; 
      std::get<1>(ctd1->CI[id]) += Wb;
      std::get<2>(ctd1->CI[id]) += delE * NN;
    }

    if (interFlag == recomb_1) {
      if (RECOMB_KE) {
	if (RECOMB_IP) {
	  lQ rQ(Z1, C1-1);
	  delE = ch.IonList[rQ]->ip + kEe1[id];
	} else
	  delE = kEe1[id];
      }

      p1->iattrib[use_key] = k1.updateC(--C1);
      partflag      = 1;
      std::get<0>(ctd1->RR[id])++; 
      std::get<1>(ctd1->RR[id]) += Wb;
      std::get<2>(ctd1->RR[id]) += delE * NN;
    }
    
    //-------------------------
    // Particle 2 interactions
    //-------------------------

    if (interFlag == free_free_2) {
      delE          = IS.selectFFInteract(ch.IonList[Q2], id);
      partflag      = 2;
      std::get<0>(ctd2->ff[id])++;
      std::get<1>(ctd2->ff[id]) += Wb;
      std::get<2>(ctd2->ff[id]) += delE * NN;
    }

    if (interFlag == colexcite_2) {
      delE         = IS.selectCEInteract(ch.IonList[Q2], CE2[id]);
      partflag     = 2;
      std::get<0>(ctd2->CE[id])++; 
      std::get<1>(ctd2->CE[id]) += Wb;
      std::get<2>(ctd2->CE[id]) += delE * NN;
    }

    if (interFlag == ionize_2) {
      delE = IS.DIInterLoss(ch.IonList[Q2]);
      p2->iattrib[use_key] = k2.updateC(++C2);
      std::get<0>(ctd2->CI[id])++; 
      std::get<1>(ctd2->CI[id]) += Wb;
      std::get<2>(ctd2->CI[id]) += delE * NN;
      partflag     = 2;
    }

    if (interFlag == recomb_2) {
      if (RECOMB_KE) {
	if (RECOMB_IP) {
	  lQ rQ(Z2, C2-1);
	  delE = ch.IonList[rQ]->ip + kEe2[id];
	} else
	  delE = kEe2[id];
      }

      p2->iattrib[use_key] = k2.updateC(--C2);
      partflag     = 2;
      std::get<0>(ctd2->RR[id])++; 
      std::get<1>(ctd2->RR[id]) += Wb;
      std::get<2>(ctd2->RR[id]) += delE * NN;
    }

    // Convert to super particle
    //
    delE *= NN;
    
    // Convert back to cgs
    //
    delE = delE * eV;
  }
  
  // For elastic interactions, delE == 0
  //
  assert(delE >= 0.0);

  // Artifically prevent cooling by setting the energy removed from
  // the COM frame to zero
  //
  if (NO_COOL) delE = 0.0;

  // Convert energy loss to system units
  //
  delE = delE/UserTreeDSMC::Eunit;

  //
  // Perform energy adjustment in ion, system COM frame with system
  // mass units
  //

  double KE1i = 0.0, KE2i = 0.0;
  double KE1f = 0.0, KE2f = 0.0;
  if (KE_DEBUG) {
    for (auto v : p1->vel) KE1i += v*v;
    for (auto v : p2->vel) KE2i += v*v;
  }

  // Mass per particle in amu for this interaction
  //
  double m1 = atomic_weights[Z1];
  double m2 = atomic_weights[Z2];

  // Total effective mass in the collision (atomic mass units)
  //
  double Mt = m1 + m2;

  // Available center of mass energy in the ballistic collision
  // (system units)
  //
  double vi  = (*cr);

  // -----------------
  // ENERGY DIAGNOSTIC
  // -----------------
  // Electrons from Particle 2 have interacted with atom/ion in Particle 1
  //
  if (partflag==1) {

    ctd1->eV_av[id] += kEe2[id];
    if (std::isnan(ctd2->eV_av[id])) {
      std::cout << "eV_N=" << ctd1->eV_N[id] << std::endl;
    }
    ctd1->eV_N[id]++;
    ctd1->eV_min[id] = std::min(ctd1->eV_min[id], kEe2[id]);
    ctd1->eV_max[id] = std::max(ctd2->eV_max[id], kEe2[id]);
    
    if (kEe2[id] > 10.2) { ctd1->eV_10[id]++;}
  }

  // -----------------
  // ENERGY DIAGNOSTIC
  // -----------------
  // Electrons from Particle 1 interacted with atom/ion in Particle 2
  //
  if (partflag==2) {

    ctd2->eV_av[id] += kEe1[id];
    if (std::isnan(ctd2->eV_av[id])) {
      std::cout << "eV_N=" << ctd2->eV_N[id] << std::endl;
    }
    ctd2->eV_N[id]++;
    ctd2->eV_min[id] = std::min(ctd2->eV_min[id], kEe1[id]);
    ctd2->eV_max[id] = std::max(ctd2->eV_max[id], kEe1[id]);
    
    if (kEe1[id] > 10.2) { ctd2->eV_10[id]++; }
  }

  // Compute the new, scattered velocities using original COM energy
  //
  std::vector<double> vrel(3), vcom(3);

  for (unsigned k=0; k<3; k++) {
    vcom[k] = (m1*p1->vel[k] + m2*p2->vel[k]) / Mt;
  }
	    
  double cos_th = 1.0 - 2.0*(*unit)();       // Cosine and sine of
  double sin_th = sqrt(1.0 - cos_th*cos_th); // Collision angle theta
  double phi    = 2.0*M_PI*(*unit)();	     // Collision angle phi
  
  vrel[0] = vi * cos_th;	  // Compute post-collision relative
  vrel[1] = vi * sin_th*cos(phi); // velocity for an elastic 
  vrel[2] = vi * sin_th*sin(phi); // interaction

				// Center of mass lost energy (m1, m2,
				// Mt are in atomic mass units)
  double deltaKE = 0.0, qKEfac = 0.5*Wa*m1*q*(1.0 - q);
  for (unsigned k=0; k<3; k++) {
    double vdif = vcom[k] + m2/Mt*vrel[k] - p1->vel[k];
    deltaKE += vdif*vdif*qKEfac;
  }

  // Update post-collision velocities.  In the electron version, the
  // momentum is assumed to be coupled to the ions, so the ion
  // momentum must be conserved.
  // 

  double kE = 0.0;		// Compute the COM KE, after scatter
  for (size_t k=0; k<3; k++) {
    p1->vel[k] = (1.0 - q)*p1->vel[k] + q*(vcom[k] + m2/Mt*vrel[k]);
    p2->vel[k] = vcom[k] - m1/Mt*vrel[k];

    double cr  = p1->vel[k] - p2->vel[k];
    kE += cr*cr;
  }

  // Reduced mass in system units for energy update
  //
  double Mu = p1->mass * p2->mass / (p1->mass + p2->mass);

  // Available KE in COM frame, system units
  //
  kE *= 0.5*Mu;

  // Total energy available in COM after removing radiative and
  // collisional loss.  A negative value for totE will be handled
  // below . . .
  //
  double totE  = kE - delE;
  double totE0 = kE - delE;

  // Cooling rate diagnostic histogram
  //
  if (TSDIAG && delE>0.0) {
				// Histogram index
    int indx = (int)floor(log(kE/delE)/(log(2.0)*TSPOW) + 5);
				// Floor and ceiling
    if (indx<0 ) indx = 0;
    if (indx>10) indx = 10;
				// Add entry
    EoverT[id][indx] += p1->mass + p2->mass;
  }
  
  //
  // Time step "cooling" diagnostic
  //
  if (use_delt>=0 && delE>0.0 && kE>0.0) {
    double dtE = kE/delE * spTau[id];
    double dt1 = p1->dattrib[use_delt];
    double dt2 = p2->dattrib[use_delt];
    p1->dattrib[use_delt] = std::max<double>(dt1, dtE);
    p2->dattrib[use_delt] = std::max<double>(dt2, dtE);
  }

  if (use_exes>=0 && delE>0.0) {
    // (-/+) value means under/overcooled: positive/negative increment
    // to delE NB: delE may be < 0 if too much energy was radiated
    // previously . . .
    //
    delE -= p1->dattrib[use_exes] + p2->dattrib[use_exes];
    p1->dattrib[use_exes] = p2->dattrib[use_exes] = 0.0;
  }
  
  lostSoFar[id] += delE;
  decelT[id]    += delE;
    
  ret            = 0;		// No error

  if (partflag==1) {
    std::get<0>(ctd1->dv[id])++; 
    std::get<1>(ctd1->dv[id]) += Wb;
    std::get<2>(ctd1->dv[id]) += delE;
  }
    
  if (partflag==2) {
    std::get<0>(ctd2->dv[id])++; 
    std::get<1>(ctd2->dv[id]) += Wb;
    std::get<2>(ctd2->dv[id]) += delE;
  }

  // Assign interaction energy variables
  //
  double Exs = 0.0;
  if (Z1 == Z2 and use_cons>=0) {
    double del = p1->dattrib[use_cons] + p2->dattrib[use_cons];
    Exs  += del;
    totE += del;
    p1->dattrib[use_cons] = p2->dattrib[use_cons] = 0.0;
  }

  // Attempt to defer negative energy adjustment
  //
  double missE = std::min<double>(0.0, totE);

  // Save energy adjustiments for next interation
  //
  if (Z1 == Z2) {		// Split between like species
    double del = 0.5*deltaKE + 0.5*missE;
    p1->dattrib[use_cons] += del;
    p2->dattrib[use_cons] += del;
  } else {				// Give excess to non-trace species
    p1->dattrib[use_cons] += deltaKE + missE;
  }

  // Compute the change of energy in the collision frame by computing
  // the velocity reduction factor
  //
  
  double vfac = 1.0;
  if (kE>0.0) vfac = totE0>0.0 ? sqrt(totE0/kE) : 0.0;

  // Update electron velocties.  Electron velocity is computed so that
  // momentum is conserved ignoring the doner ion.  Use of reduction
  // factor keeps electrons and ions in equipartition.
  //
  if (use_elec) {

    // Electron from particle #2
    //
    if (interFlag > 100 and interFlag < 200) {
      m2 = atomic_weights[0];	// Electron mass
      Mt = m1 + m2;
      for (unsigned k=0; k<3; k++) {
	vcom[k] = (m1*p1->vel[k] + m2*p2->dattrib[use_elec+k]) / Mt;
      }
	    
      vi = sqrt(2.0*kEe1[id]*eV/(atomic_weights[0]*amu)) / UserTreeDSMC::Vunit;

      vrel[0] = vi * cos_th;	      // Compute post-collision relative
      vrel[1] = vi * sin_th*cos(phi); // velocity for an elastic 
      vrel[2] = vi * sin_th*sin(phi); // interaction

      double vf2 = 0.0;		      // Debug electron energy loss/gain
      for (size_t k=0; k<3; k++) {
	double vv = p2->dattrib[use_elec+k] = vcom[k] - m1/Mt*vrel[k] * vfac;
	vf2 += vv*vv;
      }
      velER[id].push_back(vf2/(vi*vi));

    }

    // Electron from particle #1
    //
    if (interFlag > 200 and interFlag < 300) {
      m1 = atomic_weights[0];	// Electron mass
      Mt = m1 + m2;
      for(unsigned k=0; k<3; k++) {
	vcom[k] = (m1*p1->dattrib[use_elec+k] + m2*p2->vel[k]) / Mt;
      }
	    
      vi = sqrt(2.0*kEe2[id]*eV/(atomic_weights[0]*amu)) / UserTreeDSMC::Vunit;

      vrel[0] = vi * cos_th;	      // Compute post-collision relative
      vrel[1] = vi * sin_th*cos(phi); // velocity for an elastic 
      vrel[2] = vi * sin_th*sin(phi); // interaction

      double vf2 = 0.0;		      // Debug electron energy loss/gain
      for (size_t k=0; k<3; k++) {
	double vv = p1->dattrib[use_elec+k] = vcom[k] - m2/Mt*vrel[k] * vfac;
	vf2 += vv*vv;
      }
      velER[id].push_back(vf2/(vi*vi));

    }
  }

  //
  // Perform energy adjustment in ion, system COM frame with system
  // mass units.  Reduced mass has already been computed but not
  // needed below.
  //
  m1 = p1->mass;
  m2 = p2->mass;
  Mt = m1 + m2;

  for (unsigned k=0; k<3; k++) {
    vcom[k] = (m1*p1->vel[k] + m2*p2->vel[k]) / Mt;
    vrel[k] = p1->vel[k] - p2->vel[k];
  }

  vfac = 1.0;
  if (kE>0.0) vfac = totE>0.0 ? sqrt(totE/kE) : 1.0;

  *cr = 0.0;
  for (size_t k=0; k<3; k++) {
    double v1 = p1->vel[k] = vcom[k] + m2/Mt * vrel[k] * vfac;
    double v2 = p2->vel[k] = vcom[k] - m1/Mt * vrel[k] * vfac;
    double dv = p1->vel[k] - p2->vel[k];
    *cr += dv*dv;

    if (std::isnan(v1) || std::isnan(v2)) {
      std::cout << "Vel NaN" << std::endl;
    }
  }
  *cr = sqrt(*cr);

  // KE debugging
  //
  if (KE_DEBUG) {

    for (auto v : p1->vel) KE1f += v*v;
    for (auto v : p2->vel) KE2f += v*v;

				// Pre collision KE
    KE1i *= 0.5*p1->mass;
    KE2i *= 0.5*p2->mass;
				// Post collision KE
    KE1f *= 0.5*p1->mass;
    KE2f *= 0.5*p2->mass;

    double tKEi = KE1i + KE2i;	// Total pre collision KE
    double tKEf = KE1f + KE2f;	// Total post collision KE
    double dKE  = tKEi - tKEf - deltaKE; // Energy balance
    
				// Check Energy balance including excess
    if (fabs(dKE + Exs - delE - missE) > 1.0e-20)
      std::cout << "Total ("<< m1 << "," << m2 << ") = " 
		<< std::setw(14) << dKE + Exs - delE - missE
		<< ", dKE=" << std::setw(14) << dKE
		<< ", Cns=" << std::setw(14) << deltaKE
		<< ", Exs=" << std::setw(14) << Exs
		<< ", del=" << std::setw(14) << delE
		<< ", msE=" << std::setw(14) << missE
		<< ", fac=" << std::setw(14) << vfac
		<< std::endl;
  }
  
  return ret;
}


void CollideIon::debugDeltaE(double delE, unsigned short Z, unsigned short C,
			     double KE, double prob, int interFlag)
{
  if (delE < 0.0)
    std::cout << " *** Neg deltaE=" << std::setw(12) << delE << ", (Z, C)=(" 
	      << std::setw(2) << Z << ", " << std::setw(2) << C << "), E=" 
	      << std::setw(12) << KE  << ", prob=" << std::setw(12) << prob
	      << " :: " << labels[interFlag] << std::endl;
}

int CollideIon::inelasticTrace(pCell* const c, 
			       Particle* const p1, Particle* const p2, 
			       double *cr, int id)
{
  // For particle map and weights
  //
  typedef std::map<speciesKey, double> keyW;

  int ret = 0;			// No error (flag)

  // Number of protons per system mass unit
  //
  double Fn = UserTreeDSMC::Munit/amu;

  // Compute subspecies probability list
  //
  keyCross prob1, prob2;
  double   psum1 = 0.0, psum2 = 0.0;

  // Compute the fractional collision ratio for each ionic species and
  // interaction type
  //
  for (auto sp : sCross[id]) {

    for (auto v : sCross[id][sp.first]) {

      // Mass fractions
      //
      double w1 = p1->dattrib[SpList[sp.first]];
      double w2 = p2->dattrib[SpList[sp.first]];

      // Weighting by number
      //
      prob1[sp.first].push_back(v * w1 / atomic_weights[sp.first.first]);
      psum1 += v * w1 / atomic_weights[sp.first.first];

      prob2[sp.first].push_back(v * w2 / atomic_weights[sp.first.first]);
      psum2 += v * w2 / atomic_weights[sp.first.first];
    }
  }

  // Normalization for particle 1
  //
  for (auto sp : prob1) {
    for (auto &v : prob1[sp.first]) v /= psum1;
  }

  // Normalization for particle 2
  //
  for (auto sp : prob2) {
    for (auto &v : prob2[sp.first]) v /= psum2;
  }

  //-------------------------
  // VERBOSE DEBUG TEST
  //-------------------------
  // Set to false for production
  //  |
  //  v
  if (DEBUG_CR) {
    std::cout << std::setw( 8) << "index"
	      << std::setw( 4) << "Z"
	      << std::setw( 4) << "C"
	      << std::setw( 8) << "flag"
	      << std::setw(14) << "cross"
	      << std::setw(14) << "Prob 1"
	      << std::setw(14) << "Prob 2"
	      << std::setw(18) << "type label"
	      << std::endl
	      << std::setw( 8) << "-----"
	      << std::setw( 4) << "---"
	      << std::setw( 4) << "---"
	      << std::setw( 8) << "-----"
	      << std::setw(14) << "---------"
	      << std::setw(14) << "---------"
	      << std::setw(14) << "---------"
	      << std::setw(18) << "---------------"
	      << std::endl;

    // sp1 and sp2 are of type std::map<speciesKey, int>

    for (auto s : SpList) {
      speciesKey k(s.first);

      for (size_t i = 0; i < sCross[id][k].size(); i++) {
	std::cout << std::setw( 8) << i
		  << std::setw( 4) << k.first
		  << std::setw( 4) << k.second
		  << std::setw( 8) << sInter[id][k][i]
		  << std::setw(14) << sCross[id][k][i] 
		  << std::setw(14) << prob1[k][i] 
		  << std::setw(14) << prob2[k][i] 
		  << std::setw(18) << labels[sInter[id][k][i]]
		  << std::endl;
      }
    }
    std::cout << std::endl;
  }

  // Copy weights for trace components of each particle rather than
  // doing this on the fly (mostly for algorithmic clarity)
  //
  // NB: SpList maps species key to attribute position in particle
  // data
  //
  keyW new1, new2;
  for (auto sp : SpList) {
    new1[sp.first] = p1->dattrib[sp.second];
    new2[sp.first] = p2->dattrib[sp.second];
  }

  // Cycle through all pairs of species
  //
  double delE = 0.0;

  //             +--- This was computed and cached by crossSectionTrace
  //             |
  //             v
  for (auto sp : sCross[id]) {

    // The interacting species
    //
    speciesKey k0 = sp.first;

    // Number of interaction types in this map
    //
    size_t snum = sInter[id][k0].size();

    // Sanity check
    //
    if (CROSS_DBG && collD->find(k0) == collD->end()) {
      std::cout << "Missing key (" << k0.first << ", " << k0.second << ")"
		<< ", I will now crash" << std::endl 
		<< "Available keys:"  << std::endl;
      for (auto p : *collD) {
	std::cout << "** (" << p.first.first << ", " << p.first.second << ")"
		  << std::endl;
      }
      std::cout << std::endl << "Species list:" << std::endl;
      for (auto s : SpList) {
	std::cout << "** (" << s.first.first << ", " << s.first.second << ")"
		  << std::endl;
      }
    }

    // These are the collision diagnostic instances
    //
    collTDPtr ctd = (*collD)[k0];

    // Atomic number and ionic state (C=1 means neutral)
    //
    unsigned short Z = k0.first, C = k0.second;
  
    // Cycle through the interaction list
    //
    for (size_t isp=0; isp<snum; isp++) {
      
      // Get the type of this interaction
      //
      int interFlag = sInter[id][k0][isp];

      // Get trace-species mass fractions
      //
      double w1 = new1[k0];
      double w2 = new2[k0];

      // Fraction of molecules in the superparticle undergoing this
      // interaction
      //
      double Pr1 = prob1[k0][isp];
      double Pr2 = prob2[k0][isp];

      // Accumulate the total energy lost for each particle
      //
      double delE1 = 0.0, delE2 = 0.0;

      // Indicate energy loss in particle
      //
      bool pFlag = false;

      // Ion keys
      //
      lQ Q(Z, C);

      // Number of TRUE interactions in the superparticle
      //
      double N1 = Pr1 * Fn * p1->mass;
      double N2 = Pr2 * Fn * p2->mass;

      //--------------------------------------------------
      // Computation of energy and weight change for
      // each interaction type
      //--------------------------------------------------

      if (interFlag == free_free) {
	double dE = IS.selectFFInteract(ch.IonList[Q], id);

	delE1 = dE * N1;
	delE2 = dE * N2;

	std::get<0>(ctd->ff[id]) ++;
	std::get<1>(ctd->ff[id]) += Pr1 + Pr2;
	std::get<2>(ctd->ff[id]) += delE1 + delE2;

	debugDeltaE(delE1, Z, C, kEe2[id], Pr1, interFlag);
	debugDeltaE(delE2, Z, C, kEe1[id], Pr2, interFlag);

	pFlag = true;
      }

      if (interFlag == colexcite) {
	double dE = IS.selectCEInteract(ch.IonList[Q], kCE[id][k0]);

	delE1 = dE * N1;
	delE2 = dE * N2;

	std::get<0>(ctd->CE[id]) ++;
	std::get<1>(ctd->CE[id]) += Pr1 + Pr2;
	std::get<2>(ctd->CE[id]) += delE1 + delE2;

	debugDeltaE(delE1, Z, C, kCE[id][k0].back().second, Pr1, interFlag);
	debugDeltaE(delE2, Z, C, kCE[id][k0].back().second, Pr2, interFlag);

	pFlag = true;
      }

      if (interFlag == ionize) {
	double dE = IS.DIInterLoss(ch.IonList[Q]); // <--- Energy lost

	delE1 = dE * N1;
	delE2 = dE * N2;

	speciesKey kk(Z, C+1);

	double W1 = Pr1 * atomic_weights[Z];
	double W2 = Pr2 * atomic_weights[Z];

	if (W1 < w1) { // Change in mass fraction owing to ionization
	  new1[kk] += W1;
	  new1[k0] -= W1;
	} else {
	  excessW[id].push_back(dKeyD(dKey(k0, kk), p1->mass*(W1 - w1)));
	  new1[kk] += w1;
	  new1[k0]  = 0.0;
	}

	if (W2 < w2) { // Change in mass fraction owing to ionization
	  new2[kk] += W2;
	  new2[k0] -= W2;
	} else {
	  excessW[id].push_back(dKeyD(dKey(k0, kk), p2->mass*(W2 - w2)));
	  new2[kk] += w2;
	  new2[k0]  = 0.0;
	}

	std::get<0>(ctd->CI[id]) ++;
	std::get<1>(ctd->CI[id]) += Pr1 + Pr2;
	std::get<2>(ctd->CI[id]) += delE1 + delE2;

	debugDeltaE(delE1, Z, C, 0.0, Pr1, interFlag);
	debugDeltaE(delE2, Z, C, 0.0, Pr2, interFlag);

	pFlag = true;
      }

      if (interFlag == recomb) {

	if (RECOMB_KE) {
	  if (RECOMB_IP) {
	    lQ rQ(Z, C-1);
	    double Xi = ch.IonList[rQ]->ip;
	    delE1 = (Xi + kEe2[id]) * N1;
	    delE2 = (Xi + kEe1[id]) * N2;
	  } else {
	    delE1 = kEe2[id] * N1;
	    delE2 = kEe1[id] * N2;
	  }
	}

	speciesKey kk(Z, C-1);

	double W1 = Pr1 * atomic_weights[Z];
	double W2 = Pr2 * atomic_weights[Z];

	if (W1 < w1) {
	  new1[kk] += W1;
	  new1[k0] -= W1;
	} else {
	  excessW[id].push_back(dKeyD(dKey(k0, kk), p1->mass*(W1 - w1)));
	  new1[kk] += w1;
	  new1[k0]  = 0.0;
	}

	if (W2 < w2) {
	  new2[kk] += W2;
	  new2[k0] -= W2;
	} else {
	  excessW[id].push_back(dKeyD(dKey(k0, kk), p2->mass*(W2 - w2)));
	  new2[kk] += w2;
	  new2[k0]  = 0.0;
	}

	std::get<0>(ctd->RR[id]) ++;
	std::get<1>(ctd->RR[id]) += Pr1 + Pr2;
	std::get<2>(ctd->RR[id]) += delE1 + delE2;

	debugDeltaE(delE1, Z, C, kEe2[id], Pr1, interFlag);
	debugDeltaE(delE2, Z, C, kEe1[id], Pr2, interFlag);

	pFlag = true;
      }

      // Energy diagnostics
      //
      if (pFlag) {
	double kEe = kEe1[id];

	ctd->eV_av[id] += kEe;
	if (std::isnan(ctd->eV_av[id])) {
	  std::cout << "eV_N=" << ctd->eV_N[id] << std::endl;
	}
	ctd->eV_N[id]++;
	ctd->eV_min[id] = std::min(ctd->eV_min[id], kEe);
	ctd->eV_max[id] = std::max(ctd->eV_max[id], kEe);
	
	if (kEe > 10.2) { ctd->eV_10[id]++; }

	if (N1>0.0) {
	  std::get<0>(ctd->dv[id]) ++; 
	  std::get<1>(ctd->dv[id]) += Pr1;
	  std::get<2>(ctd->dv[id]) += delE1/N1;
	}

	if (N2>0.0) {
	  std::get<0>(ctd->dv[id]) ++; 
	  std::get<1>(ctd->dv[id]) += Pr2;
	  std::get<2>(ctd->dv[id]) += delE2/N2;
	}
      }

      // Convert back to cgs
      //
      delE += (delE1 + delE2) * eV;

      if (delE < 0.0) {
	std::cout << "Found delE=" << std::setw(14) << delE/UserTreeDSMC::Eunit 
		  << ", delE1=" << std::setw(14) << delE1
		  << ", delE2=" << std::setw(14) << delE2
		  << ", w1="    << std::setw(14) << w1
		  << ", w2="    << std::setw(14) << w2
		  << ", N1="    << std::setw(14) << N1
		  << ", N2="    << std::setw(14) << N2
		  << ", (Z, C) = (" 
		  << std::setw(2) << Z << ", "
		  << std::setw(2) << C << ") "
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
  double kE = 0.5*Mu*(*cr)*(*cr);

  // For tracking energy conservation (system units)
  //
  double dE = kE*TolV*TolV;
    
  // Artifically prevent cooling by setting the energy removed from
  // the COM frame to zero
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
    double dtE = remE/delE * spTau[id];
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
  
  double norm1 = 0.0, norm2 = 0.0;
  for (auto sp : SpList) {
    norm1 += new1[sp.first];
    norm2 += new2[sp.first];
  }

  if (fabs(norm1-1.0)>1.0e-8) {
    std::cout << "Norm1 error: " << norm1-1.0 << std::endl;
    for (auto sp : SpList) new1[sp.first] /= norm1;
  }

  if (fabs(norm2-1.0)>1.0e-8) {
    std::cout << "Norm2 error: " << norm2-1.0 << std::endl;
    for (auto sp : SpList) new2[sp.first] /= norm2;
  }

  // Replace particle weights
  //
  for (auto sp : SpList) {
    p1->dattrib[sp.second] = new1[sp.first];
    p2->dattrib[sp.second] = new2[sp.first];
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
  *epsm    = 0.0;		// EPSM to be implemented . . . 
}


void * CollideIon::timestep_thread(void * arg)
{
  pHOT* tree = (pHOT* )((tstep_pass_arguments*)arg)->tree;
  int id     = (int)((tstep_pass_arguments*)arg)->id;
  
  thread_timing_beg(id);
  
  // Loop over cells, cell time-of-flight time for each particle
  //
  pCell *c;
  Particle *p;
  double L, DT, mscale;
  
  // Loop over cells, processing collisions in each cell
  //
  for (unsigned j=0; j<cellist[id].size(); j++ ) {
    
    // Number of particles in this cell
    //
    c = cellist[id][j];
    L = c->Scale();
 
    double volc = c->Volume();
    
    sKeyDmap            densM, lambdaM, crossM;
    sKey2Dmap           crossIJ;
    
    crossIJ = totalScatteringCrossSections(0, c, id);
    
    for (auto it1 : c->count) {
      speciesKey i1 = it1.first;
      densM[i1] = c->Mass(i1)/volc;
    }
    
    double meanDens=0.0, meanLambda=0.0;
    
    if (MFPTS) {

      for (auto it1 : c->count) {
	speciesKey i1 = it1.first;
	crossM [i1] = 0.0;

	for (auto it2 : c->count) {
	  
	  speciesKey i2 = it2.first;
	  double      N = UserTreeDSMC::Munit/amu;
	  
	  if (i2 == defaultKey) N /= atomic_weights[i2.first];
	  
	  if (i2>=i1) {
	    crossM[i1] += N * densM[i2] * crossIJ[i1][i2];
	  } else
	    crossM[i1] += N * densM[i2] * crossIJ[i2][i1];
	}
	
	lambdaM[i1] = 1.0/crossM[i1];
	meanDens   += densM[i1] ;
	meanLambda += densM[i1] * lambdaM[i1];
      }
    
      // This is the number density-weighted
      meanLambda /= meanDens;
    }
    
    for (auto i : c->bods) {

      // Current particle
      //
      p = tree->Body(i);

      // Compute time of flight criterion and assign cell scale to
      // characteristic size
      //
      DT     = 1.0e40;
      mscale = 1.0e40;

      double vtot = 0.0;
      for (unsigned k=0; k<3; k++) {
	mscale = std::min<double>(pHOT::sides[k]*L, mscale);
	vtot += p->vel[k]*p->vel[k];
      }
      vtot = sqrt(vtot) + 1.0e-40;
      
      // Compute collision time criterion
      //
      if (MFPTS) {
	for (unsigned k=0; k<3; k++)
	  DT = std::min<double>(meanLambda/vtot, DT);
      }

      // Size scale for multistep timestep calc.
      //
      p->scale = mscale;

      // Compute cooling criterion timestep
      //
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

void CollideIon::finalize_cell(pHOT* const tree, pCell* const cell, 
			       double kedsp, int id)
{
  //
  // Spread out species change differences
  //
  if (aType == Trace) {

    if (excessW[id].size()) {

      std::vector<unsigned long> bods = cell->bods;

      //
      // Sanity check [prior]
      //
      if (EXCESS_DBG) {
	for (auto b : bods) {
	  Particle *p = tree->Body(b);
	  double sum  = 0.0;
	  for (auto s : SpList) {
				// Check value and accumulate
	    if (std::isnan(p->dattrib[s.second]))
	      std::cout << "[" << myid
			<< "] body #" << b
			<< ": NaN weight! (prior)" << std::endl;
	    if (std::isinf(p->dattrib[s.second]))
	      std::cout << "[" << myid
			<< "] body #" << b
			<< ": inf weight! (prior)" << std::endl;
	    if (p->dattrib[s.second]<0.0) {
	      std::cout << "[" << myid
			<< "] body #" << b
			<< ": negative weight! (prior)" << std::endl;
	    } else {
	      sum += p->dattrib[s.second];
	    }
	  }
				// Check final sum
	  if (fabs(sum - 1.0) < 1.0e-12) {
	    std::cout << "[" << myid
		      << "] body #" << b
		      << ": normalization error=" << sum << " (prior)" << std::endl;
	  }
	}

	for (auto w : excessW[id]) {
	  if (std::isnan(w.second)) {
	    speciesKey k1 = w.first.first;
	    speciesKey k2 = w.first.second;
	    std::cout << "[" << myid
		      << "] excess NaN "
		      << "(" << k1.first << ", " << k1.second << ")"
		      << "(" << k2.first << ", " << k2.second << ")"
		      << std::endl;
	  }
	}

      } // end: Sanity check


      keyWghts totW, sumW;
      for (auto s : SpList) totW[s.first] = sumW[s.first] = 0.0;
				// Get the total cell mass in each trace sp
      for (auto b : bods) {
	Particle *p = tree->Body(b);
	for (auto s : SpList) 
	  totW[s.first] += p->mass * p->dattrib[s.second];
      }

				// Get the transfer excess in each trace sp
      for (auto w : excessW[id]) sumW[w.first.first] += w.second;

      if (EXCESS_DBG) {
	for (auto s : SpList) {
	  if (std::isnan(totW[s.first])) {
	    std::cout << "[" << myid
		      << "] totW is NaN "
		      << "(" << s.first.first << ", " << s.first.second << ")"
		      << std::endl;
	  }
	  if (std::isinf(totW[s.first])) {
	    std::cout << "[" << myid
		      << "] totW is Inf "
		      << "(" << s.first.first << ", " << s.first.second << ")"
		      << std::endl;
	  }
	  if (std::isnan(sumW[s.first])) {
	    std::cout << "[" << myid
		      << "] sumW is NaN "
		      << "(" << s.first.first << ", " << s.first.second 
		      << ")" << std::endl;
	    for (auto w : excessW[id]) {
	      speciesKey k1 = w.first.first;
	      speciesKey k2 = w.first.second;
	      std::cout << "(" << k1.first << "," << k1.second << "):"
			<< "(" << k2.first << "," << k2.second << ")"
			<< "-->" << w.second << std::endl;
	    }
	  }
	  if (std::isinf(sumW[s.first])) {
	    std::cout << "[" << myid
		      << "] sumW is Inf "
		      << "(" << s.first.first << ", " << s.first.second 
		      << ")" << std::endl;
	    for (auto w : excessW[id]) {
	      speciesKey k1 = w.first.first;
	      speciesKey k2 = w.first.second;
	      std::cout << "(" << k1.first << "," << k1.second << "):"
			<< "(" << k2.first << "," << k2.second << ")"
			<< "-->" << w.second << std::endl;
	    }
	  }
	}
      }

				// Process the excess list
      for (auto w : excessW[id]) {
				// Remove FROM this species
	speciesKey k1 = w.first.first;
				// Add TO this species
	speciesKey k2 = w.first.second;
				// More removals than mass?
	w.second *= std::min<double>(1.0, totW[k1]/sumW[k1]);
				// Shuffle the body indices
	std::random_shuffle(bods.begin(), bods.end());
				// Loop through the particles
	for (auto b : bods) {
	  Particle *p = tree->Body(b);
	
	  int j1 = SpList[k1];	// <From> species index
	  int j2 = SpList[k2];	// <To  > species index

	  // Skip if particle doesn't have this trace species
	  //
	  if (p->dattrib[j1] > 0.0) {
	    double ww = w.second/p->mass;
	    if (EXCESS_DBG) {
	      if (std::isnan(ww)) {
		std::cout << "[" << myid
			  << "] weight is NaN, mass=" << p->mass
			  << "(" << k1.first << ", " << k1.second << ")"
			  << "(" << k2.first << ", " << k2.second << ")"
			  << std::endl;
	      }
	      if (std::isinf(ww)) {
		std::cout << "[" << myid
			  << "] weight is Inf, mass=" << p->mass
			  << "(" << k1.first << ", " << k1.second << ")"
			  << "(" << k2.first << ", " << k2.second << ")"
			  << std::endl;
	      }
	    }

	    if (ww > p->dattrib[j1]) {
	      w.second -= p->mass * p->dattrib[j1];
	      p->dattrib[j2] += p->dattrib[j1];
	      p->dattrib[j1]  = 0.0;
	    } else {
	      w.second = 0.0;
	      p->dattrib[j2] += ww;
	      p->dattrib[j1] -= ww;
	    }
	  }
	  
	  if (w.second<=0) break;
	}
      }

      //
      // Sanity check
      //
      if (EXCESS_DBG) {
	for (auto b : bods) {
	  Particle *p = tree->Body(b);
	  double sum  = 0.0;
	  for (auto s : SpList) {
				// Check value and accumulate
	    if (std::isnan(p->dattrib[s.second]))
	      std::cout << "[" << myid
			<< "] body #" << b
			<< ": NaN weight! (posterior)" << std::endl;
	    if (std::isinf(p->dattrib[s.second]))
	      std::cout << "[" << myid
			<< "] body #" << b
			<< ": inf weight! (posterior)" << std::endl;
	    if (p->dattrib[s.second]<0.0) {
	      std::cout << "[" << myid
			<< "] body #" << b
			<< ": negative weight! (posterior)" << std::endl;
	    } else {
	      sum += p->dattrib[s.second];
	    }
	  }
				// Check final sum
	  if (fabs(sum - 1.0) < 1.0e-12) {
	    std::cout << "[" << myid
		      << "] body #" << b
		      << ": normalization error=" << sum 
		      << " (posterior)" << std::endl;
	  }
	}

      } // end: Sanity check

    }
  }
  
  //
  // Cross-section debugging
  //
  if (CROSS_DBG && id==0) {
    if (nextTime_dbg <= tnow && nCnt_dbg < nCel_dbg)
      {
	speciesKey i;

	if (aType==Direct or aType==Weight)
	  i = cell->count.begin()->first;
	else
	  i = SpList.begin()->first;

	cross2_dbg.push_back(csections[id][i][i]);
	nCnt_dbg++;
	if (nCnt_dbg == nCel_dbg) write_cross_debug();
      }
  }

  //
  // Done
  //
}

// Help class that maintains database of diagnostics
//
collDiag::collDiag(CollideIon* caller) : p(caller)
{
  // Initialize the map
  //
  if (p->ZList.size()) {

    for (auto n : p->ZList) {

      unsigned short Z = n;

      for (unsigned short C=1; C<Z+2; C++) {
	speciesKey k(Z, C);
	(*this)[k] = collTDPtr(new CollisionTypeDiag());
      }
    }
  } else if (p->SpList.size()) {
    for (auto n : p->SpList) {
      (*this)[n.first] = collTDPtr(new CollisionTypeDiag());
    }
  } else {			// Sanity check
    if (myid==0) {
      std::cerr << "collDiag:CollDiag: species list or map is "
		<< "not initialized" << std::endl;
    }
    MPI_Abort(MPI_COMM_WORLD, 57);
  }

  Esum.resize(nthrds, 0.0);

  // Initialize the output file
  //
  initialize();
}

// Gather statistics from all processes
//
void collDiag::gather()
{
  for (auto it : *this) {
    collTDPtr ctd = it.second;
    ctd->sumUp();
    ctd->sync();

    Esum_s = std::accumulate(Esum.begin(), Esum.end(), 0.0);
    double z;
    MPI_Reduce(&(z=Esum_s), &Esum_s, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

// Zero out the counters
//
void collDiag::reset()
{
  for (auto it : *this) it.second->reset();
  std::fill(Esum.begin(), Esum.end(), 0.0);
}

void collDiag::initialize() 
{
  if (myid) return;

  {
    // Generate the file name
    std::ostringstream sout;
    sout << outdir << runtag << ".ION_coll";
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
	    << "# W(ff)         summed wght of free-free " << std::endl
	    << "# E(ff)         cum energy in free-free  " << std::endl
	    << "# N(ce)         number of collisions     " << std::endl
	    << "# W(ce)         summed wght of collision " << std::endl
	    << "# E(ce)         cum energy in collisions " << std::endl
	    << "# N(ci)         number of ionizations    " << std::endl
	    << "# W(ci)         summed wght of ionized   " << std::endl
	    << "# E(ci)         cum energy in ionizations" << std::endl
	    << "# N(rr)         number of rad recombs    " << std::endl
	    << "# W(rr)         summed wght of recombs   " << std::endl
	    << "# E(rr)         energy in rad recombs    " << std::endl
	    << "# d(KE)         mean energy change       " << std::endl
	    << "# Elost         total energy loss        " << std::endl
	    << "# Etotl         total kinetic energy     " << std::endl
	    << "#"                                         << std::endl;
	
				// Species labels
	out << "#" << std::setw(11+12) << std::right << "Species==>" << " | ";
	for (auto it : *this) {
	  ostringstream sout, sout2;
	  sout  << "(" << it.first.first << ", " << it.first.second << ")";
	  size_t w =13*12, l = sout.str().size();
	  sout2 << std::setw((w-l)/2) << ' ' << sout.str();
	  out   << std::setw(w) << sout2.str() << " | ";
	}
	out << std::setw(2*12) << ' ' << " |" << std::endl;

				// Header line
	out << std::setfill('-') << std::right;
	out << "#" << std::setw(11+12) << '+' << " | ";
	for (auto it : *this) {
	  for (int i=0; i<13; i++) out << std::setw(12) << '+';
	  out << " | ";
	}
	out << std::setw(12) << '+' << std::setw(12) << '-' << " |"
	    << std::setfill(' ') << std::endl;

				// Column labels
	out << "#" 
	    << std::setw(11) << "Time |"
	    << std::setw(12) << "Temp |" << " | ";
	for (auto it : *this) {
	  out << std::setw(12) << "N(ff) |"
	      << std::setw(12) << "W(ff) |"
	      << std::setw(12) << "E(ff) |"
	      << std::setw(12) << "N(ce) |"
	      << std::setw(12) << "W(ce) |"
	      << std::setw(12) << "E(ce) |"
	      << std::setw(12) << "N(ci) |"
	      << std::setw(12) << "W(ci) |"
	      << std::setw(12) << "E(ci) |"
	      << std::setw(12) << "N(rr) |"
	      << std::setw(12) << "W(rr) |"
	      << std::setw(12) << "E(rr) |"
	      << std::setw(12) << "d(KE) |"
	      << " | ";
	}
	out << std::setw(12) << "Elost |" << std::setw(12) << "Etotl"
	    << " |" << std::endl;
	
				// Column numbers
	std::ostringstream st;
	unsigned int cnt = 0;
	st << "[" << ++cnt << "] |";
	out << "#" << std::setw(11) << st.str();
	st.str("");
	st << "[" << ++cnt << "] |";
	out << std::setw(12) << st.str() << " | ";
	for (auto it : *this) {
	  for (size_t l=0; l<13; l++) {
	    st.str("");
	    st << "[" << ++cnt << "] |";
	    out << std::setw(12) << std::right << st.str();
	  }
	  out << " | ";
	}
	st.str("");
	st << "[" << ++cnt << "] |";
	out << std::setw(12) << std::right << st.str();
	st.str("");
	st << "[" << ++cnt << "]";
	out << std::setw(12) << std::right << st.str() << " |" << std::endl;

				// Header line
	out << std::setfill('-') << std::right;
	out << "#" << std::setw(11+12) << '+' << " | ";
	for (auto it : *this) {
	  for (int i=0; i<13; i++) out << std::setw(12) << '+';
	  out << " | ";
	}
	out << std::setw(12) << '+' << std::setw(12) << '-' << " |"
	    << std::setfill(' ') << std::endl;
      }
    }
    in.close();
  }

  {
    // Generate the file name
    std::ostringstream sout;
    sout << outdir << runtag << ".ION_energy";
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
	for (auto it : *this) {
	  ostringstream sout, sout2;
	  sout  << "(" << it.first.first << ", " << it.first.second << ")";
	  size_t w = 5*12, l = sout.str().size();
	  sout2 << std::setw((w-l)/2) << ' ' << sout.str();
	  out   << std::setw(w) << sout2.str() << " | ";
	}
	out << std::endl;
	
				// Header line
	out << std::setfill('-') << std::right;
	out << "#" << std::setw(11) << '+' << " | ";
	for (auto it : *this) {
	  for (int i=0; i<5; i++) out << std::setw(12) << '+';
	  out << " | ";
	}
	out << std::setfill(' ') << std::endl;

				// Column labels
	out << "#" << std::setw(11) << "Time" << " | ";
	for (auto it : *this) {
	  out << std::setw(12) << "avg |"
	      << std::setw(12) << "num |"
	      << std::setw(12) << "min |"
	      << std::setw(12) << "max |"
	      << std::setw(12) << "over10 |"
	      << " | ";
	}
	out << std::endl;

				// Column numbers
	std::ostringstream st;
	unsigned int cnt = 0;
	st << "[" << ++cnt << "] |";
	out << "#" << std::setw(11) << st.str() << " | ";
	for (auto it : *this) {
	  for (size_t l=0; l<5; l++) {
	    st.str("");
	    st << "[" << ++cnt << "] |";
	    out << std::setw(12) << std::right << st.str();
	  }
	  out << " | ";
	}
	out << std::endl;
				// Header line
	out << std::setfill('-') << std::right;
	out << "#" << std::setw(11) << '+' << " | ";
	for (auto it : *this) {
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
      double Etot   = 0.0;
      double cvrt   = eV/UserTreeDSMC::Eunit;

      out << std::setw(12) << tnow 
	  << std::setw(12) << p->tempM << " | ";
      for (auto it : *this) {
	collTDPtr ctd = it.second;
	out << std::setw(12) << std::get<0>(ctd->ff_s)
	    << std::setw(12) << std::get<1>(ctd->ff_s)
	    << std::setw(12) << std::get<2>(ctd->ff_s) * cvrt
	    << std::setw(12) << std::get<0>(ctd->CE_s)
	    << std::setw(12) << std::get<1>(ctd->CE_s)
	    << std::setw(12) << std::get<2>(ctd->CE_s) * cvrt
	    << std::setw(12) << std::get<0>(ctd->CI_s)
	    << std::setw(12) << std::get<1>(ctd->CI_s)
	    << std::setw(12) << std::get<2>(ctd->CI_s) * cvrt
	    << std::setw(12) << std::get<0>(ctd->RR_s)
	    << std::setw(12) << std::get<1>(ctd->RR_s)
	    << std::setw(12) << std::get<2>(ctd->RR_s) * cvrt;
	if (std::isnan(std::get<2>(ctd->dv_s))) {
	  std::cout << "Error in print, t=" << tnow << std::endl;
	}
	if (std::get<1>(ctd->dv_s)>0.0)
	  out << std::setw(12) << std::get<2>(ctd->dv_s)/std::get<1>(ctd->dv_s) << " | ";
	else
	  out << std::setw(12) << 0.0 << " | ";
	Etot += 
	  std::get<2>(ctd->ff_s) + std::get<2>(ctd->CE_s) +
	  std::get<2>(ctd->CI_s) + std::get<2>(ctd->RR_s) ;
      }
      out << std::setw(12) << Etot * cvrt 
	  << std::setw(12) << Esum_s
	  << " |" << std::endl;
    }
  }

  {
    std::ofstream out(energy_file_debug.c_str(), ios::out | ios::app);
    if (out) {
      out << std::setw(12) << tnow << " | ";
      for (auto it : *this) {
	collTDPtr ctd = it.second;
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
  unsigned char nOK = 0;
  
  //
  // Default values
  //
  use_cons = -1;
  use_elec = -1;

  //
  // Let root node ONLY do the reading
  //
  if (myid==0) {

    std::ifstream in(map.c_str());
    if (in.bad()) 
      {
	std::cerr << "CollideIon::parseSpecies: definition file <"
		  << map << "> could not be opened . . . quitting"
		  << std::endl;
	nOK = 1;
      }


    // Ok, read the first line to get the implementation type

    if (nOK == 0) {

      const int nline = 2048;
      char line[nline];

      in.getline(line, nline);
      std::string type(line);

      if (type.compare("direct")==0) {
    
	aType = Direct;

	if (use_key<0) {
	  std::cerr << "CollideIon: species key position is not defined in "
		    << "Component" << std::endl;
	  nOK = 1;
	}

	if (nOK == 0) {
	  
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
	}

      } else if (type.compare("weight")==0) {
    
	aType = Weight;

	if (use_key<0) {
	  std::cerr << "CollideIon: species key position is not defined in "
		    << "Component" << std::endl;
	  nOK = 1;
	}

	in.getline(line, nline);
	
	if (in.good()) {
	  std::istringstream sz(line);
	  sz >> use_cons;
	  if (sz.good()) {
	    sz >> use_elec;
	  }
	} else {
	  nOK = 1;		// Can't read electrons or use_cons value, fatal
	}

				// Print warning, not fatal
	if (use_cons<0) {
	  std::cout << "CollideIon: energy key position is not defined, "
		    << "you using trace-species weighting but not imposing energy conservation"
		    << std::endl;
	}

	if (use_elec<0) {
	  std::cout << "CollideIon: electron key position is not defined, "
		    << "you using trace-species weighting with electron velocity emulation"
		    << std::endl;
	}

	if (nOK == 0) {
	  
	  int Z;
	  double W, M;
	  while (1) {
	    in.getline(line, nline);
	    if (in.good()) {
	      std::istringstream sz(line);
	      sz >> Z;
	      sz >> W;
	      sz >> M;		// Add to the element list
	      if (!sz.bad()) {
		ZList.insert(Z);
		ZWList[Z] = W;
		ZMList[Z] = M;
	      }
	    } else {
	      break;
	    }
	  }

	  // Find the largest weight (assume fiducial)
	  double sMax = 0.0;
	  for (auto v : ZWList) {
	    if (v.second > sMax) {
	      sFid = v.first;
	      sMax = v.second;
	    }
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
	std::cerr << "CollideIon::parseSpecies: implementation type <"
		  << type << "> is not recognized . . . quitting"
		  << std::endl;
	nOK = 1;
      }
    }

  }

  MPI_Bcast(&nOK, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

  
  if (nOK) MPI_Abort(MPI_COMM_WORLD, 55);

  int is = aType;
  MPI_Bcast(&is, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (myid) {
    switch (is) {
    case Direct:
      aType = Direct;
      break;
    case Weight:
      aType = Weight;
      break;
    case Trace:
      aType = Trace;
      break;
    default:
      std::cout << "Proc " << myid << ": error in enum <" << is << ">"
		<< std::endl;
      MPI_Abort(MPI_COMM_WORLD, 56);
    }
  }
      
  unsigned int sz = ZList.size();
  MPI_Bcast(&sz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  unsigned short z;

  if (myid==0) {
    for (auto it : ZList) {
      z = it;
      MPI_Bcast(&z, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
    }
  } else {
    for (unsigned j=0; j<sz; j++) {
      MPI_Bcast(&z, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
      ZList.insert(z);
    }
  }


  sz = ZWList.size();
  MPI_Bcast(&sz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  double v;

  if (myid==0) {
    for (auto it : ZWList) {
      z = it.first;
      v = it.second;
      MPI_Bcast(&z, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&v, 1, MPI_DOUBLE,         0, MPI_COMM_WORLD);
    }

    for (auto it : ZMList) {
      z = it.first;
      v = it.second;
      MPI_Bcast(&z, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&v, 1, MPI_DOUBLE,         0, MPI_COMM_WORLD);
    }

  } else {
    for (unsigned j=0; j<sz; j++) {
      MPI_Bcast(&z, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&v, 1, MPI_DOUBLE,         0, MPI_COMM_WORLD);
      ZWList[z] = v;
    }

    for (unsigned j=0; j<sz; j++) {
      MPI_Bcast(&z, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&v, 1, MPI_DOUBLE,         0, MPI_COMM_WORLD);
      ZMList[z] = v;
    }
  }

  MPI_Bcast(&use_cons, 1, MPI_INT,         0, MPI_COMM_WORLD);
  MPI_Bcast(&use_elec, 1, MPI_INT,         0, MPI_COMM_WORLD);

  if (aType == Trace) {

    speciesKey key;
    int pos;

    sz = SpList.size();
    MPI_Bcast(&sz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (myid==0) {
      for (auto it : SpList) {

	key = it.first;
	pos = it.second;

	MPI_Bcast(&key.first,  1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&key.second, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&pos,        1, MPI_INT,            0, MPI_COMM_WORLD);
      }
    } else {
      for (unsigned j=0; j<sz; j++) {
	MPI_Bcast(&key.first,  1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&key.second, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&pos,        1, MPI_INT,            0, MPI_COMM_WORLD);
	SpList[key] = pos;
      }
    }
  }
  
}

sKey2Umap CollideIon::generateSelection
(pCell* const c, sKeyDmap* const Fn, double crm, double tau, int id,
 double& meanLambda, double& meanCollP, double& totalNsel)
{
  if (aType == Direct)
    return generateSelectionDirect(c, Fn, crm, tau, id, 
				   meanLambda, meanCollP, totalNsel);
  else if (aType == Weight)
    return generateSelectionWeight(c, Fn, crm, tau, id, 
				   meanLambda, meanCollP, totalNsel);
  else
    return generateSelectionTrace(c, Fn, crm, tau, id, 
				  meanLambda, meanCollP, totalNsel);
}

sKey2Umap CollideIon::generateSelectionDirect
(pCell* const c, sKeyDmap* const Fn, double crm, double tau, int id,
 double& meanLambda, double& meanCollP, double& totalNsel)
{
  sKeyDmap            densM, collPM, lambdaM, crossM;
  sKey2Dmap           selcM;
  sKey2Umap           nselM;
    
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
  //
  // Done
  //
  
  for (auto it1 : c->count) {
    speciesKey i1 = it1.first;
    densM[i1] = c->Mass(i1)/volc / atomic_weights[i1.first];
    //                             ^
    //                             |
    // Number density--------------+
    //
  }
    
  double meanDens = 0.0;
  meanLambda      = 0.0;
  meanCollP       = 0.0;
    
  for (auto it1 : c->count) {

    speciesKey i1 = it1.first;
    crossM [i1]   = 0.0;

    for (auto it2 : c->count) {

      speciesKey i2 = it2.first;

      if (i2>=i1) {
	crossM[i1] += (*Fn)[i2]*densM[i2]*csections[id][i1][i2];
      } else
	crossM[i1] += (*Fn)[i2]*densM[i2]*csections[id][i2][i1];
      
      if (csections[id][i1][i2] <= 0.0 || std::isnan(csections[id][i1][i2])) {
	cout << "INVALID CROSS SECTION! :: " << csections[id][i1][i2]
	     << " #1 = (" << i1.first << ", " << i1.second << ")"
	     << " #2 = (" << i2.first << ", " << i2.second << ")";
	csections[id][i1][i2] = 0.0; // Zero out
      }
	    
      if (csections[id][i2][i1] <= 0.0 || std::isnan(csections[id][i2][i1])) {
	cout << "INVALID CROSS SECTION! :: " << csections[id][i2][i1]
	     << " #1 = (" << i2.first << ", " << i2.second << ")"
	     << " #2 = (" << i1.first << ", " << i1.second << ")";
	csections[id][i2][i1] = 0.0; // Zero out
      }
	
    }
      
    if (it1.second>0 && (crossM[i1] == 0 || std::isnan(crossM[i1]))) {
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
    
  // This is the number density-weighted MFP (used for diagnostics
  // only)
  //
  meanLambda /= meanDens;

  // Number-density weighted collision probability (used for
  // diagnostics only)
  //
  meanCollP  /= meanDens;
    
  // This is the per-species N_{coll}
  //
  totalNsel = 0.0;

  std::map<speciesKey, unsigned>::iterator it1, it2;

  for (it1=c->count.begin(); it1!=c->count.end(); it1++) {
    speciesKey i1 = it1->first;
    
    for (it2=it1; it2!=c->count.end(); it2++) {
      speciesKey i2 = it2->first;
      
      // Probability of an interaction of between particles of type 1
      // and 2 for a given particle of type 2
      //
      double Prob = (*Fn)[i2] * densM[i2] * csections[id][i1][i2] * crm * tau;
      
      // Count _pairs_ of identical particles only
      //                 |
      //                 v
      if (i1==i2)
	selcM[i1][i2] = 0.5 * (it1->second-1) *  Prob;
      else
	selcM[i1][i2] = it1->second * Prob;
      //
      // For double-summing of species A,B and B,A interactions 
      // when A != B is list orders A<B and therefore does not double 
      // count (see line 951 in Collide.cc)
      
      nselM[i1][i2] = static_cast<unsigned>(floor(selcM[i1][i2]+0.5));
      totalNsel += nselM[i1][i2];
    }
  }
  
  return nselM;
}

sKey2Umap CollideIon::generateSelectionWeight
(pCell* const c, sKeyDmap* const Fn, double crm, double tau, int id,
 double& meanLambda, double& meanCollP, double& totalNsel)
{
  sKeyDmap            densM, densN, collP, nsigmaM, ncrossM;
  sKey2Dmap           selcM;
  sKey2Umap           nselM;
    
  // Convert from CHIANTI to system units
  //
  const double cunit = 1e-14/(UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);

  // Sample cell
  //
  pCell *samp = c->sample;

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
  //
  // Done
  //
  
  {
    for (auto it1 : c->count) {

      // Only compute if particles of this species is in the cell
      //
      if (it1.second) {
	speciesKey i1 = it1.first;

	// Mass density scaled by atomic weight in amu
	//
	densM[i1] = c->Mass(i1) / atomic_weights[i1.first] / volc;
	
	// Number density of superparticles
	//
	densN[i1] = static_cast<double>(c->Count(i1))/volc;
      }
    }
  }
    
  if (0) {
    std::cout << std::endl
	      << std::setw(10) << "Species"
	      << std::setw(16) << "n dens"
	      << std::setw(16) << "m dens"
	      << std::setw(16) << "sp mass"
	      << std::setw(10) << "n count"
	      << std::endl
	      << std::setw(10) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(10) << "---------"
	      << std::endl;

    for (auto it : c->count) {
      std::ostringstream sout;
      sout << "(" << it.first.first << ", " << it.first.second << ")";
      std::cout << std::setw(10) << sout.str()
		<< std::setw(16) << densN[it.first]
		<< std::setw(16) << densM[it.first]
		<< std::setw(16) << c->Mass(it.first)
		<< std::setw(10) << c->bods.size()
		<< std::endl;
    }
  }

  if (DEBUG_SL) {
    std::cout << std::endl
	      << std::setw(16) << "Species"
	      << std::setw(16) << "Cross"
	      << std::setw(16) << "elec V"
	      << std::setw(16) << "densM"
	      << std::setw(16) << "densN"
	      << std::setw(10) << "count 1"
	      << std::setw(10) << "count 2"
	      << std::endl
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(10) << "---------"
	      << std::setw(10) << "---------"
	      << std::endl;

    std::map<speciesKey, unsigned>::iterator it1, it2;
    
    for (it1=c->count.begin(); it1!=c->count.end(); it1++) {

      if (it1->second==0) continue;
      
      speciesKey k1 = it1->first;

      for (it2=it1; it2!=c->count.end(); it2++) {

	if (it2->second==0) continue;

	speciesKey k2 = it2->first;

	std::ostringstream sout;
	sout << "<" 
	     << k1.first << "," << k1.second << "|"
	     << k2.first << "," << k2.second << ">";
	std::cout << std::setw(16) << sout.str()
		  << std::setw(16) << csections[id][k1][k2] / cunit
		  << std::setw(16) << Evel[id]
		  << std::setw(16) << densM[k1]
		  << std::setw(16) << densN[k1]
		  << std::setw(10) << it1->second
		  << std::setw(10) << it2->second
		  << std::endl;
      }
    }
  }

  double meanDens = 0.0;
  meanLambda      = 0.0;
  meanCollP       = 0.0;
    
  for (auto it1 : c->count) {

    // Only compute if particles of this species is in the cell
    //
    if (it1.second) {

      speciesKey i1 = it1.first;
      ncrossM[i1]   = 0.0;
      nsigmaM[i1]   = 0.0;

      for (auto it2 : c->count) {

	// Only compute if particles of this species is in the cell
	//
	if (it2.second) {

	  speciesKey i2 = it2.first;

	  // Compute the computational cross section (that is, true
	  // cross seciton scaled by number of true particles per
	  // computational particle)

	  double crossT = 0.0;
	  if (i2>=i1)
	    crossT = csections[id][i1][i2];
	  else
	    crossT = csections[id][i2][i1];

	  // Choose the trace species of the two (may be neither in
	  // which case it doesn't matter)

	  if (densM[i2] <= densM[i1]) {
	    unsigned Z2  = i2.first;
	    crossT      *= (*Fn)[i2] * ZMList[Z2] / atomic_weights[Z2];
	    ncrossM[i1] += crossT;
	    nsigmaM[i1] += densN[i2]*crossT;
	  } else {
	    unsigned Z1  = i1.first;
	    crossT      *= (*Fn)[i1] * ZMList[Z1] / atomic_weights[Z1];
	    ncrossM[i2] += crossT;
	    nsigmaM[i2] += densN[i1]*crossT;
	  }
      
	  // So, ncrossM is the superparticle cross section for each species

	  if (csections[id][i1][i2] <= 0.0 || std::isnan(csections[id][i1][i2])) {
	    cout << "INVALID CROSS SECTION! :: " << csections[id][i1][i2]
		 << " #1 = (" << i1.first << ", " << i1.second << ")"
		 << " #2 = (" << i2.first << ", " << i2.second << ")";

	    csections[id][i1][i2] = 0.0; // Zero out
	  }
	  
	  if (csections[id][i2][i1] <= 0.0 || std::isnan(csections[id][i2][i1])) {
	    cout << "INVALID CROSS SECTION! :: " << csections[id][i2][i1]
		 << " #1 = (" << i2.first << ", " << i2.second << ")"
		 << " #2 = (" << i1.first << ", " << i1.second << ")";

	    csections[id][i2][i1] = 0.0; // Zero out
	  }
	}
      }
    }
  }
      
  for (auto it1 : c->count) {

    // Only compute if particles of this species is in the cell
    //
    if (it1.second) {

      speciesKey i1 = it1.first;

      if (ncrossM[i1] == 0 || std::isnan(ncrossM[i1])) {
	cout << "INVALID CROSS SECTION! ::"
	     << " (" << i1.first << ", " << i1.second << ")"
	     << " nsigmaM = " << nsigmaM [i1]
	     << " ncrossM = " << ncrossM [i1] 
	     << " Fn = "      <<   (*Fn) [i1] << endl;
      
	std::cout << std::endl
		  << std::setw(10) << "Species"
		  << std::setw(16) << "x-section"
		  << std::setw(16) << "sp mass"
		  << std::setw(16) << "n*sigma"
		  << std::setw(16) << "n*cross"
		  << std::endl
		  << std::setw(10) << "---------"
		  << std::setw(16) << "---------"
		  << std::setw(16) << "---------"
		  << std::setw(16) << "---------"
		  << std::setw(16) << "---------"
		  << std::endl;

	for (auto it : csections[id][i1]) {
	  std::ostringstream sout;
	  sout << "(" << it.first.first << ", " << it.first.second << ")";
	  cout << std::setw(10) << sout.str()
	       << std::setw(16) << it.second 
	       << std::setw(16) << c->Mass(it.first)
	       << std::setw(16) << nsigmaM[it.first]
	       << std::setw(16) << ncrossM[it.first] 
	       << std::endl;
	}
      }
    
      collP  [i1] = nsigmaM[i1] * crm * tau;
    
      meanDens   += densN[i1];
      meanCollP  += densN[i1] * collP  [i1];
      meanLambda += densN[i1] * nsigmaM[i1];
    }
  }
    

  if (DEBUG_SL) {
    std::cout << std::endl
	      << std::setw(10) << "Species"
	      << std::setw(16) << "count"
	      << std::setw(16) << "sp mass"
	      << std::setw(16) << "n*sigma"
	      << std::setw(16) << "n*cross"
	      << std::setw(16) << "Prob"
	      << std::endl
	      << std::setw(10) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::endl;

    for (auto it : c->count) {

      // Only output if particles of this species is in the cell
      //
      if (it.second) {
	double prob = densN[it.first] * ncrossM[it.first] * crm * tau;
	std::ostringstream sout;
	sout << "(" << it.first.first << ", " << it.first.second << ")";
	cout << std::setw(10) << sout.str()
	     << std::setw(16) << it.second 
	     << std::setw(16) << c->Mass(it.first)
	     << std::setw(16) << nsigmaM[it.first]
	     << std::setw(16) << ncrossM[it.first]
	     << std::setw(16) << prob
	     << std::endl;
      }
    }
  }

  // This is the number density-weighted MFP (used for diagnostics
  // only)
  //
  meanLambda  = meanDens/meanLambda;

  // Number-density weighted collision probability (used for
  // diagnostics only)
  //
  meanCollP  /= meanDens;
    
  // This is the per-species N_{coll}
  //
  totalNsel = 0.0;

  std::map<speciesKey, unsigned>::iterator it1, it2;

  for (it1=c->count.begin(); it1!=c->count.end(); it1++) {

    // Only compute if particles of this species is in the cell
    if (it1->second) {

      speciesKey i1 = it1->first;
    
      for (it2=it1; it2!=c->count.end(); it2++) {

	// Only compute if particles of this species is in the cell
	if (it2->second) {

	  speciesKey i2 = it2->first;
      
	  double crsvel = 0.0;

	  sKeyPair k(i1, i2);
	  if (i1>=i2) k = sKeyPair(i2, i1);

	  if (samp)
	    crsvel = std::get<0>(ntcdb[samp->mykey]->VelCrsAvg(k, 0.95));
	  else
	    crsvel = std::get<0>(ntcdb[c->mykey]->VelCrsAvg(k, 0.95));
	  
	  // Probability of an interaction of between particles of type 1
	  // and 2 for a given particle of type 2
	  //
	  double Prob = 0.0;

	  if (densM[i1]>=densM[i2]) {
	    Prob = densM[i2] * (*Fn)[i2] * cunit * crsvel * tau;
	  } else {
	    Prob = densM[i1] * (*Fn)[i1] * cunit * crsvel * tau;
	  }

	  // Count _pairs_ of identical particles only
	  //                 |
	  //                 |
	  if (i1==i2) //     v
	    selcM[i1][i2] = 0.5 * it1->second * (it2->second-1) *  Prob;
	  else
	    selcM[i1][i2] = it1->second * it2->second * Prob;
	
	  // For debugging only
	  //
	  if (DEBUG_SL) {
	    if (selcM[i1][i2]>10000) {
	      double cv1, cv2, cv3;
	      if (samp) {
		cv1 = std::get<0>(ntcdb[samp->mykey]->VelCrsAvg(k, 0.50));
		cv2 = std::get<0>(ntcdb[samp->mykey]->VelCrsAvg(k, 0.90));
		cv3 = std::get<0>(ntcdb[samp->mykey]->VelCrsAvg(k, 0.95));
	      } else {
		cv1 = std::get<0>(ntcdb[c->mykey]->VelCrsAvg(k, 0.50));
		cv2 = std::get<0>(ntcdb[c->mykey]->VelCrsAvg(k, 0.90));
		cv3 = std::get<0>(ntcdb[c->mykey]->VelCrsAvg(k, 0.95));
	      }

	      std::cout << "Too many collisions: collP=" << meanCollP
			<< ", MFP=" << meanLambda << ", P=" << Prob
			<< ", <sigma*vel>=" << crsvel
			<< ", N=" << selcM[i1][i2]
			<< ", q(0.5, 0.9, 0.95) = (" << cv1 << ", "
			<< cv2 << ", " << cv3 << ")"
			<< std::endl;
	    }
	  }

	  //
	  // For double-summing of species A,B and B,A interactions 
	  // when A != B is list orders A<B and therefore does not double 
	  // count (see line 951 in Collide.cc)
	  
	  nselM[i1][i2] = static_cast<unsigned>(floor(selcM[i1][i2]+0.5));
	  totalNsel += nselM[i1][i2];
	}
      }
    }
  }
  

  if (0) {
    unsigned nbods  = c->bods.size();
    double totalEst = 0.5 * meanCollP * nbods * (nbods-1);
    if (totalNsel > 200.0 && totalNsel > 5.0 * totalEst) {
      std::cout << "Total pairs: " << totalNsel  << std::endl
		<< "  Est pairs: " << totalEst   << std::endl
		<< "     mean P: " << meanCollP  << std::endl
		<< "     bodies: " << nbods      << std::endl;
    }
  }

  if (DEBUG_SL) {

    std::cout << std::endl
	      << std::endl     << std::right
	      << std::setw(16) << "Interact"
	      << std::setw(16) << "N sel"
	      << std::setw(16) << "Prob 0"
	      << std::setw(16) << "Prob 1"
	      << std::endl
	      << std::setw(16) << "--------"
	      << std::setw(16) << "--------"
	      << std::setw(16) << "--------"
	      << std::setw(16) << "--------"
	      << std::endl;
      
    for (it1=c->count.begin(); it1!=c->count.end(); it1++) {

      // Only output if particles of this species is in the cell
      //
      if (it1->second) {

	for (it2=it1; it2!=c->count.end(); it2++) {
	  
	  // Only output if particles of this species is in the cell
	  //
	  if (it2->second) {

	    speciesKey i1 = it1->first;
	    speciesKey i2 = it2->first;
	    sKeyPair   k(i1, i2);
      
	    double crsvel = 0.0;
	    if (samp)
	      crsvel = std::get<0>(ntcdb[samp->mykey]->VelCrsAvg(k, 0.95));
	    else
	      crsvel = std::get<0>(ntcdb[c->mykey]->VelCrsAvg(k, 0.95));

	    double Prob0 = 0.0, Prob1 = 0.0;

	    if (densM[i1]>=densM[i2]) {
	      Prob0 = densM[i2] * (*Fn)[i2] * cunit * crsvel * tau;
	      Prob1 = nsigmaM[i2] * crm * tau;
	    } else {
	      Prob0 = densM[i1] * (*Fn)[i1] * cunit * crsvel * tau;
	      Prob1 = nsigmaM[i1] * crm * tau;
	    }
	    
	    std::cout << "(" 
		      << std::setw(2) << i1.first << ","
		      << std::setw(2) << i1.second << ") ("
		      << std::setw(2) << i2.first << ","
		      << std::setw(2) << i2.second << ")  "
		      << std::setw(16) << selcM[i1][i2]
		      << std::setw(16) << Prob0
		      << std::setw(16) << Prob1
		      << std::endl;
	  }
	}
      }
    }
    std::cout << std::endl 
	      << "  Mean Coll P = " << meanCollP 
	      << "  Mean Lambda = " << meanLambda
	      << "  totalNsel = "   << totalNsel
	      << std::endl << std::endl;
  }
  
  return nselM;
}

sKey2Umap CollideIon::generateSelectionTrace
(pCell* const c, sKeyDmap* const Fn, double crm, double tau, int id,
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
  //
  if (std::isnan(csections[id][key][key]) or csections[id][key][key] < 0.0) {
    cout << "[" << myid << "] INVALID CROSS SECTION! :: " 
	 << csections[id][key][key] << std::endl;
    
    // Verbose debugging
    //
    if (EXCESS_DBG) {
      cout << "[" << myid << "] SpList size :: " << SpList .size() << std::endl;
      cout << "[" << myid << "] C body size :: " << c->bods.size() << std::endl;

      keyWghts mW;
      for (auto s : SpList) mW[s.first] = 0.0;
      double massP = 0.0;
      for (auto b : c->bods) {
	Particle *p = curTree->Body(b);
	massP += p->mass;
	for (auto s : SpList) 
	  mW[s.first] += p->mass * p->dattrib[s.second];
      }
      for (auto s : SpList)
	std::cout << std::setw(3) << s.first.first 
		  << std::setw(3) << s.first.second
		  << " : "
		  << std::setw(18) << mW[s.first]/massP
		  << std::endl;
    }

    // Zero out
    //
    csections[id][key][key] = 0.0;
  }
    
  // Cache relative velocity
  //
  spCrm[id] = crm;

  // Compute collision rates in system units
  //
  double crossM = (*Fn)[key] * dens * csections[id][key][key];
  double collPM = crossM * crm * tau;

  // Interaction rate
  //
  double rateF = (*Fn)[key] * crm * tau;

  // Cache probability of an interaction of between the particles pair
  // for use in inelasticTrace
  //
  spProb[id] = dens * rateF * num * 
    1e-14 / (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);

  // Cache time step for estimating "over" cooling timestep is use_delt>=0
  //
  spTau[id]  = tau;

  // For Collide diagnostics
  //
  meanLambda = 1.0/crossM;
  meanCollP  = collPM;
    
  double Prob  = dens * rateF * csections[id][key][key];
  double selcM = 0.5 * (num-1) * Prob;
  //             ^      ^
  //             |      |
  //             |      +--- For correct Poisson statistics
  //             |
  //             +--- Pairs are double counted
  //

  sKey2Umap nselM;
  nselM[key][key] = static_cast<unsigned>(floor(selcM+0.5));
  spNsel[id] = nselM[key][key];
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
  cross1_dbg.clear();
  cross2_dbg.clear();
}



void CollideIon::gatherSpecies()
{
  const double Tfac = 2.0*UserTreeDSMC::Eunit/3.0 * amu  /
    UserTreeDSMC::Munit/boltz;

  if (aType==Direct or aType==Weight) {

    // Compute temperature only

    double mass  = 0.0;

    consE = 0.0;
    totlE = 0.0;
    tempM = 0.0;
    tempE = 0.0;

    typedef std::map<key_type, double> cType;
    cType ETcache;

    // Interate through all cells
    //
    pHOT_iterator itree(*c0->Tree());
    
    while (itree.nextCell()) {
      
      pCell *cell = itree.Cell();
      
      // Compute the mass-weighted temerature and mass
      //
      double KEtot, KEdsp;
      cell->sample->KE(KEtot, KEdsp);

      double T = KEdsp * Tfac * molWeight(cell->sample);
      
      mass  += cell->Mass();
      tempM += cell->Mass() * T;
      totlE += cell->Mass() * KEtot;

      if (std::isnan(KEtot)) {
	std::cout << "NaN" << std::endl;
	cell->sample->KE(KEtot, KEdsp);
      }

      if (aType==Weight and use_cons >= 0) {

	for (auto b : cell->bods) {
	  consE += c0->Tree()->Body(b)->dattrib[use_cons];
	}
	
	cType::iterator ft = ETcache.find(cell->sample->mykey);

	if (ft != ETcache.end()) {
	  T = ft->second;
	} else {
	  std::vector<double> vel1(3, 0.0), vel2(3, 0.0);
	  double count = 0.0;
	  for (auto c : cell->sample->children) {
	    for (auto b : c.second->bods) {
	      Particle *p = c0->Tree()->Body(b);
	      KeyConvert k(p->iattrib[use_key]);
	      unsigned short ne = k.C() - 1;
	      double numb = p->mass/atomic_weights[k.Z()] * ne;
	      if (ne) {
		for (unsigned k=0; k<3; k++) {
		  double v = p->dattrib[use_elec+k];
		  vel1[k] += v   * numb;
		  vel2[k] += v*v * numb;
		}
		count   += numb;
	      }
	    }
	  }

	  double dispr = 0.0;
	  if (count > 0.0) {
	    for (unsigned k=0; k<3; k++) 
	      dispr += 0.5*(vel2[k] - vel1[k]*vel1[k]/count);
	    T = ETcache[cell->sample->mykey] = 
	      dispr/count * Tfac * atomic_weights[0];
	  } else {
	    T = ETcache[cell->sample->mykey] = 0.0;
	  }
	}
	tempE += cell->Mass() * T;
      }
    }


    // Send values to root
    //
    double val1, val2, val3 = 0.0, val4 = 0.0, val5 = 0.0;
    
    for (int i=1; i<numprocs; i++) {

      if (i == myid) {
				// Mass
	MPI_Send(&mass,  1, MPI_DOUBLE, 0, 331, MPI_COMM_WORLD);
				// Temp
	MPI_Send(&tempM, 1, MPI_DOUBLE, 0, 332, MPI_COMM_WORLD);

	if (aType==Weight and use_cons >= 0) {
	  MPI_Send(&consE, 1, MPI_DOUBLE, 0, 333, MPI_COMM_WORLD);
	  MPI_Send(&totlE, 1, MPI_DOUBLE, 0, 334, MPI_COMM_WORLD);
	  if (use_elec >= 0)
	    MPI_Send(&tempE, 1, MPI_DOUBLE, 0, 335, MPI_COMM_WORLD);
	}

      }
				// Root receives from Node i
      if (0 == myid) {

	MPI_Recv(&val1, 1, MPI_DOUBLE, i, 331, MPI_COMM_WORLD, 
		 MPI_STATUS_IGNORE);
	MPI_Recv(&val2, 1, MPI_DOUBLE, i, 332, MPI_COMM_WORLD, 
		 MPI_STATUS_IGNORE);

	if (aType==Weight and use_cons >= 0) {
	  MPI_Recv(&val3, 1, MPI_DOUBLE, i, 333, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  MPI_Recv(&val4, 1, MPI_DOUBLE, i, 334, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  if (use_elec >= 0)
	    MPI_Recv(&val5, 1, MPI_DOUBLE, i, 335, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);
	}

	if (std::isnan(val3) || std::isnan(val4) || std::isnan(val5)) {
	  std::cout << "NaN" << std::endl;
	}
	  

	mass  += val1;
	tempM += val2;
	consE += val3;
	totlE += val4;
	tempE += val5;
      }
    }

    if (mass>0.0) {
      tempM /= mass;
      tempE /= mass;
    }
    
  } else {

    // Clean the maps
    //
    double mass = 0.0;
    
    tempM = 0.0;
    specM.clear();

    // Interate through all cells
    //
    pHOT_iterator itree(*c0->Tree());
    
    while (itree.nextCell()) {
      
      pCell *cell = itree.Cell();
      
      // Compute the temerature
      //
      double KEtot, KEdsp;
      cell->sample->KE(KEtot, KEdsp);

      double T = KEdsp * Tfac * molWeight(cell->sample);
      
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

    // Send the temperature and local map to root
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

	MPI_Recv(&sizm, 1, MPI_INT,    i, 330, MPI_COMM_WORLD, 
		 MPI_STATUS_IGNORE);
	MPI_Recv(&val1, 1, MPI_DOUBLE, i, 331, MPI_COMM_WORLD, 
		 MPI_STATUS_IGNORE);
	MPI_Recv(&val2, 1, MPI_DOUBLE, i, 332, MPI_COMM_WORLD, 
		 MPI_STATUS_IGNORE);

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
}
  

// Print out species counts
//
void CollideIon::printSpecies
(std::map<speciesKey, unsigned long>& spec, double T)
{
  if (myid) return;

  if (aType == Direct) {	// Call the generic printSpecies member
    Collide::printSpecies(spec, tempM);
  } else if (aType == Weight) {	// Call the weighted printSpecies version
    printSpeciesWeight(spec, tempM);
  } else {			// Call the trace fraction version
    printSpeciesTrace();
  }

}

void CollideIon::electronGather()
{
  if (aType==Weight && use_elec >= 0) {

    std::vector<double> eVel, iVel;

    // Interate through all cells
    //
    pHOT_iterator itree(*c0->Tree());
    
    while (itree.nextCell()) {
      
      for (auto b : itree.Cell()->bods) {
	double cri = 0.0, cre = 0.0;
	for (int l=0; l<3; l++) {
	  double ve = c0->Tree()->Body(b)->dattrib[use_elec+l];
	  cre += ve*ve;
	  double vi = c0->Tree()->Body(b)->vel[l];
	  cri += vi*vi;
	}
	eVel.push_back(sqrt(cre));
	iVel.push_back(sqrt(cri));
      }
    }

    // Accumulate from threads
    //
    std::vector<double> loss;
    for (int t=0; t<nthrds; t++) {
      loss.insert(loss.end(), velER[t].begin(), velER[t].end());
      velER[t].clear();
    }

    for (int i=1; i<numprocs; i++) {

      if (i == myid) {
	unsigned eNum = eVel.size();
	MPI_Send(&eNum,       1, MPI_UNSIGNED, 0, 335, MPI_COMM_WORLD);
	MPI_Send(&eVel[0], eNum, MPI_DOUBLE,   0, 336, MPI_COMM_WORLD);
	MPI_Send(&iVel[0], eNum, MPI_DOUBLE,   0, 337, MPI_COMM_WORLD);

	eNum = loss.size();
	MPI_Send(&eNum,       1, MPI_UNSIGNED, 0, 338, MPI_COMM_WORLD);
	MPI_Send(&loss[0], eNum, MPI_DOUBLE,   0, 339, MPI_COMM_WORLD);
      }
				// Root receives from Node i
      if (0 == myid) {
	unsigned eNum;
	MPI_Recv(&eNum,       1, MPI_UNSIGNED, i, 335, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	std::vector<double> vTmp(eNum);
	MPI_Recv(&vTmp[0], eNum, MPI_DOUBLE,   i, 336, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	eVel.insert(eVel.begin(), vTmp.begin(), vTmp.end());
	MPI_Recv(&vTmp[0], eNum, MPI_DOUBLE,   i, 337, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	iVel.insert(iVel.begin(), vTmp.begin(), vTmp.end());
	MPI_Recv(&eNum,       1, MPI_UNSIGNED, i, 338, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	vTmp.resize(eNum);
	MPI_Recv(&vTmp[0], eNum, MPI_DOUBLE,   i, 339, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	loss.insert(loss.end(), vTmp.begin(), vTmp.end());
      }
    }
    
    if (myid==0) {

      if (eVel.size()) {
	// Sort the lists
	std::sort(eVel.begin(), eVel.end());
	std::sort(iVel.begin(), iVel.end());

	// Make the histograms
	elecH = ahistoDPtr(new AsciiHisto<double>(eVel, 20, 0.01));
	ionH  = ahistoDPtr(new AsciiHisto<double>(iVel, 20, 0.01));

	// Make the quantiles
	size_t qnt_s = qnt.size(), ev_s = eVel.size();
	elecV.resize(qnt_s);
	ionV .resize(qnt_s);
	for (size_t i=0; i<qnt_s; i++) {
	  elecV[i] = eVel[floor(ev_s*qnt[i])];
	  ionV [i] = iVel[floor(ev_s*qnt[i])];
	}
      }

      if (loss.size()) {
	lossH = ahistoDPtr(new AsciiHisto<double>(loss, 20, 0.01));
      }
    }

  }
}
  

void CollideIon::electronPrint(std::ostream& out)
{
  // Print the header for electron quantiles
  //
  out << std::endl << std::string(53, '-')  << std::endl
      << "-----Electron velocity quantiles---------------------" << std::endl
      << std::string(53, '-') << std::endl << std::left
      << std::setw(12) << "Quantile" 
      << std::setw(16) << "V_electron"
      << std::setw(16) << "V_ion"      << std::endl
      << std::setw(12) << "--------" 
      << std::setw(16) << "----------"
      << std::setw(16) << "----------" << std::endl;
  for (size_t i=0; i<qnt.size(); i++)
    out << std::setw(12) << qnt[i] 
	<< std::setw(16) << elecV[i]
	<< std::setw(16) << ionV [i] << std::endl;
  if (elecH.get()) {
    out << std::endl
	<< std::string(53, '-')  << std::endl
	<< "-----Electron velocity distribution------------------" << std::endl
	<< std::string(53, '-')  << std::endl;
    (*elecH)(out);
  }
  if (ionH.get()) {
    out << std::endl
	<< std::string(53, '-')  << std::endl
	<< "-----Ion velocity distribution-----------------------" << std::endl
	<< std::string(53, '-')  << std::endl;
    (*ionH)(out);
  }
  if (lossH.get()) {
    out << std::endl
	<< std::string(53, '-')  << std::endl
	<< "-----Electron energy gain/loss distribution----------" << std::endl
	<< std::string(53, '-')  << std::endl;
    (*lossH)(out);
  }
  out << std::endl << std::endl;
}

const std::string clabl(unsigned c)
{
  std::ostringstream sout;
  sout << "[" << c << "]  ";
  return sout.str();
}

// Print out species counts (Trace version)
//
void CollideIon::printSpeciesTrace()
{
  std::ofstream dout;

  // Generate the file name, if it does not exist
  //
  if (species_file_debug.size()==0) {
    std::ostringstream sout;
    sout << outdir << runtag << ".species";
    species_file_debug = sout.str();

    // Check for existence of file
    //
    std::ifstream in (species_file_debug.c_str());

    // Write a new file?
    //
    if (in.fail()) {

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
    }
  }

  // Open for append
  //
  if (!dout.is_open()) 
    dout.open(species_file_debug.c_str(), ios::out | ios::app);

  dout << std::setprecision(5);
  dout << "  " 
       << std::setw(12) << std::right << tnow
       << std::setw(12) << std::right << tempM;
  for (spDItr it=specM.begin(); it != specM.end(); it++)
    dout << std::setw(12) << std::right << it->second;
  dout << std::endl;
}


// Compute the mean molecular weight in atomic mass units
//
double CollideIon::molWeight(sCell *cell)
{
  double mol_weight = 1.0;

  if (aType==Direct or aType==Weight) {
    double numbC = 0.0, massC = 0.0;
    for (auto it : cell->count) {
      speciesKey i = it.first;
      double M = cell->Mass(i);
      numbC += M / atomic_weights[i.first];
      massC += M;
    }

    mol_weight = massC/numbC;
  }

  if (0) {
    double numbC = 0.0, massC = 0.0;
    for (auto it : cell->count) {
      speciesKey i = it.first;
      double M = cell->Mass(i);
      numbC += M * ZWList[i.first];
      massC += M;
    }

    mol_weight = massC/numbC;
  }

  return mol_weight;
}


void CollideIon::printSpeciesWeight(std::map<speciesKey, unsigned long>& spec,
				    double temp)
{
  if (myid) return;

  typedef std::map<speciesKey, unsigned long> spCountMap;
  typedef spCountMap::iterator spCountMapItr;

				// Field width
  const unsigned short wid = 16;

  std::ofstream dout;

				// Generate the file name
  if (species_file_debug.size()==0) {
    std::ostringstream sout;
    sout << outdir << runtag << ".species";
    species_file_debug = sout.str();

    // Check for existence of file
    //
    std::ifstream in (species_file_debug.c_str());
    
    // Write a new file?
    //
    if (in.fail()) {
      
      // Open the file for the first time
      //
      dout.open(species_file_debug.c_str());

      // Print the header
      //
      dout << "# " 
	   << std::setw(wid) << std::right << "Time "
	   << std::setw(wid) << std::right << "Temp ";
      for (spCountMapItr it=spec.begin(); it != spec.end(); it++) {
	std::ostringstream sout;
	sout << "(" << it->first.first << "," << it->first.second << ") ";
	dout << setw(wid) << right << sout.str();
      }
      if (use_cons>=0) {
	dout << std::setw(wid) << std::right << "Cons_E"
	     << std::setw(wid) << std::right << "Totl_E"
	     << std::setw(wid) << std::right << "Comb_E";
	if (use_elec>=0)
	  dout << std::setw(wid) << std::right << "Temp_E";
      }
      dout << std::endl;
      
      dout << "# " 
	   << std::setw(wid) << std::right << "--------"
	   << std::setw(wid) << std::right << "--------";
      for (spCountMapItr it=spec.begin(); it != spec.end(); it++)
	dout << setw(wid) << std::right << "--------";
      if (use_cons>=0) {
	dout << std::setw(wid) << std::right << "--------"
	     << std::setw(wid) << std::right << "--------"
	     << std::setw(wid) << std::right << "--------";
	if (use_elec>=0)
	  dout << std::setw(wid) << std::right << "--------";
      }
      dout << std::endl;
      
    }
  }

  // Open for append
  //
  if (!dout.is_open())
    dout.open(species_file_debug.c_str(), ios::out | ios::app);


  double tmass = 0.0;
  for (spCountMapItr it=spec.begin(); it != spec.end(); it++)
    tmass += ZMList[it->first.first] * it->second;

				// Use total mass to print mass
				// fraction
  dout << "  " 
       << std::setw(wid) << std::right << tnow
       << std::setw(wid) << std::right << temp;

  for (spCountMapItr it=spec.begin(); it != spec.end(); it++) {
    if (tmass > 0.0) 
      dout << std::setw(wid) << std::right 
	   << ZMList[it->first.first] * it->second / tmass;
    else
      dout << std::setw(wid) << std::right << 0.0;
  }
  if (use_cons>=0) {
    dout << std::setw(wid) << std::right << consE
	 << std::setw(wid) << std::right << totlE
	 << std::setw(wid) << std::right << totlE + consE;
    if (use_elec>=0)
      dout << std::setw(wid) << std::right << tempE;
  }
  dout << std::endl;
}

