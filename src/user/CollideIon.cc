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

double   CollideIon::Nmin          = 1.0e-08;	   
double   CollideIon::Nmax    	   = 1.0e+25;	   
double   CollideIon::Tmin    	   = 1.0e+03;	   
double   CollideIon::Tmax    	   = 1.0e+08;	   
double   CollideIon::TolV    	   = 1.0e-03;	   
unsigned CollideIon::Nnum    	   = 400;	   
unsigned CollideIon::Tnum    	   = 200;	   
string   CollideIon::cache   	   = ".HeatCool";
bool     CollideIon::frost_warning = false; // For debugging . . . 

// Artifically prevent cooling by setting the energy removed from the
// COM frame to zero
//
bool NO_COOL           = false;

// Subtract KE from COM pair for testing only.  This is technically
// incorrect since the electrons are "trace" species and not part of
// the energy conservation.
//
const bool RECOMB_KE   = true;

// Cross-section debugging; set to false for production
//
const bool CROSS_DBG   = false;

// Excess trace map debugging; set to false for production
//
const bool EXCESS_DBG  = false;


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
  spTau    .resize(nthrds);
  spCrm    .resize(nthrds);
  spNsel   .resize(nthrds);
  spProb   .resize(nthrds);

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
void CollideIon::initialize_cell(pHOT* tree, pCell* cell, double rvmax, int id)
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

  // Representative avg cell energy in ergs
  //
  double Eerg = 0.5*vavg*vavg*amu;

  // In eV
  //
  double EeV = Eerg / eV;

  // Upscaling factor for scattering cross section
  //
  double sUp  = diamfac*diamfac;

  if (aType == Direct or aType == Weight) {

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

	double CrossG = 0.5*M_PI*(Radius1 + Radius2)*(Radius1 + Radius2);
	double Cross1 = CrossG, Cross2 = CrossG;

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
    // double eVel = sqrt(amu*meanM[id]/me);
    double eVel = 1.0;

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
	Cross += M_PI*b*b * eVel * meanE[id];
      }

      double tCross = Cross * sUp * 1e-14 / 
	(UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);

      csections[id][defaultKey][defaultKey] += tCross * 
	meanF[id][k] /atomic_weights[k.first];
    }
  }
}


sKey2Dmap& 
CollideIon::totalScatteringCrossSections(double crm, pCell *c, int id)
{
  // it1 and it2 are of type std::map<speciesKey, unsigned>
  
  double vel  = crm * UserTreeDSMC::Vunit;
  double Eerg = 0.5*vel*vel*amu;
  double EeV  = Eerg / eV;

  // Upscaling factor for scattering cross section
  //
  double sUp  = diamfac*diamfac;
  
  if (aType == Direct or aType == Weight) {

    for (auto it1 : c->count) {

      speciesKey i1 = it1.first;
      double geom1  = geometric(i1.first);
    
      for (auto it2 : c->count) {

	speciesKey i2 = it2.first;
	double geom2  = geometric(i2.first);
	  
	double mu = atomic_weights[i1.first] * atomic_weights[i1.first] / 
	  (atomic_weights[i1.first] + atomic_weights[i2.first]);
	  
	double eVel1 = sqrt(amu*atomic_weights[i1.first]/me);
	double eVel2 = sqrt(amu*atomic_weights[i2.first]/me);
	  
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
    // In auto iterators below:
    //    s is of type std::map<speciesKey, int>
    //    b is of type std::vector<unsigned long>
    //
				// Equipartition electron velocity
				//
    // double eVel = sqrt(amu*meanM[id]/me);
    double eVel = 1.0;

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
	Cross += M_PI*b*b * eVel * meanE[id];
      }
      
      double tCross = Cross * sUp * 1e-14 / 
	(UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
	
      csections[id][defaultKey][defaultKey] += tCross * 
	meanF[id][k]/atomic_weights[k.first];
    }
  }

  return csections[id];
}

double CollideIon::crossSectionDirect(pHOT *tree, Particle* p1, Particle* p2, 
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
  
  // Save the per-interaction cross sections
  dCross[id].clear();

  // Index the interactions
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
    dCross[id].push_back(cross12*diamfac*diamfac);
    dInter[id].push_back(neut_neut_1);

    cross21 = geometric(Z2);
    dCross[id].push_back(cross21*diamfac*diamfac);
    dInter[id].push_back(neut_neut_2);
  }

				//-------------------------------
				// Electrons in second particle
				//-------------------------------
  if (ne2 > 0) {
    if (C1==1) {		// Neutral atom-electron scattering
      cross12 = elastic(Z1, kEe2[id]) * eVel2*ne2;
      dCross[id].push_back(cross12);
      dInter[id].push_back(neut_elec_1);
    }  else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C1-1) /
	std::max<double>(kEe2[id]*eV, FloorEv*eV) * 1.0e7; // nm
      cross12 = M_PI*b*b * eVel2*ne2;
      dCross[id].push_back(cross12);
      dInter[id].push_back(ion_elec_1);
    }
  }
    
				//-------------------------------
				// Electrons in first particle
				//-------------------------------
  if (ne1 > 0) {
    if (C2==1) {		// Neutral atom-electron scattering
      cross21 = elastic(Z2, kEe1[id]) * eVel1*ne1;
      dCross[id].push_back(cross21);
      dInter[id].push_back(neut_elec_2);
    } else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C2-1) /
	std::max<double>(kEe1[id]*eV, FloorEv*eV) * 1.0e7; // nm
      cross21 = M_PI*b*b * eVel1*ne1;
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
    double crs = eVel2*ne2 * ff1;

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

    double crs = eVel2*ne2 * CE1[id].back().first;

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
    double crs = eVel2*ne2 * DI1;

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
    double crs = eVel2*ne2 * RE1.back();

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
    double crs = eVel1*ne1 * ff2;

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
    double crs = eVel1*ne1 * CE2[id].back().first;

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
    double crs = eVel1*ne1*RE2.back();

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

double CollideIon::crossSectionWeight(pHOT *tree, Particle* p1, Particle* p2, 
				      double cr, int id)
{
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
  double N1 = p1->mass*UserTreeDSMC::Munit/amu * ZWList[Z1];
  double N2 = p2->mass*UserTreeDSMC::Munit/amu * ZWList[Z2];

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
  
  // Save the per-interaction cross sections
  dCross[id].clear();

  // Index the interactions
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
    dCross[id].push_back(cross12*diamfac*diamfac);
    dInter[id].push_back(neut_neut_1);

    cross21 = geometric(Z2);
    dCross[id].push_back(cross21*diamfac*diamfac);
    dInter[id].push_back(neut_neut_2);
  }

				//-------------------------------
				// Electrons in second particle
				//-------------------------------
  if (ne2 > 0) {
    if (C1==1) {		// Neutral atom-electron scattering
      cross12 = elastic(Z1, kEe2[id]) * eVel2*ne2;
      dCross[id].push_back(cross12);
      dInter[id].push_back(neut_elec_1);
    }  else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C1-1) /
	std::max<double>(kEe2[id]*eV, FloorEv*eV) * 1.0e7; // nm
      cross12 = M_PI*b*b * eVel2*ne2;
      dCross[id].push_back(cross12);
      dInter[id].push_back(ion_elec_1);
    }
  }
    
				//-------------------------------
				// Electrons in first particle
				//-------------------------------
  if (ne1 > 0) {
    if (C2==1) {		// Neutral atom-electron scattering
      cross21 = elastic(Z2, kEe1[id]) * eVel1*ne1;
      dCross[id].push_back(cross21);
      dInter[id].push_back(neut_elec_2);
    } else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C2-1) /
	std::max<double>(kEe1[id]*eV, FloorEv*eV) * 1.0e7; // nm
      cross21 = M_PI*b*b * eVel1*ne1;
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
    double crs = eVel2*ne2 * ff1;

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

    double crs = eVel2*ne2 * CE1[id].back().first;

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
    double crs = eVel2*ne2 * DI1;

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
    double crs = eVel2*ne2 * RE1.back();

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
    double crs = eVel1*ne1 * ff2;

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
    double crs = eVel1*ne1 * CE2[id].back().first;

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
    double crs = eVel1*ne1*RE2.back();

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


double CollideIon::crossSectionTrace(pHOT *tree, Particle* p1, Particle* p2, 
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
  // double eVel = sqrt(meanM[id]*amu/me);
  double eVel = 1.0;

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

				// Upscaling factor (default: 1)
      double sUp = diamfac * diamfac;
				// Geometric cross sections based on
				// atomic radius
      double Radius = geometric(Z) + meanR[id];

      crossS = neutF[id] * M_PI*Radius*Radius * sUp;
      tCrossMap.push_back(crossS);
      tInterMap.push_back(neut_neut);
	
      //
      // Neutral atom-electron scattering
      //
      
      crossS = elastic(Z, kEe) * eVel * meanE[id];
      
      tCrossMap.push_back(crossS);
      tInterMap.push_back(neut_elec);
      
    } else {
				// 
				// Rutherford scattering
				//
      double b = 0.5*esu*esu*(C-1) /
	std::max<double>(kEe*eV, FloorEv*eV) * 1.0e7; // nm
      crossS = M_PI*b*b * eVel * meanE[id];

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


int CollideIon::inelasticDirect(pHOT *tree, Particle* p1, Particle* p2, 
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
  double dE   = kE*TolV*TolV;
  double delE  = 0.0, delEeV = 0.0;

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
    //          |
    //          v
    const bool DEBUG_F = false;
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
	lQ rQ(Z1, C1-1);
	delE = ch.IonList[rQ]->ip;
      }

      p1->iattrib[use_key] = k1.updateC(--C1);
      partflag      = 1;
      std::get<0>(ctd1->RR[id])++; 
      std::get<1>(ctd1->RR[id]) += NN;
      std::get<2>(ctd1->RR[id]) += kEe2[id] * NN;
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
	lQ rQ(Z1, C1-1);
	delE = ch.IonList[rQ]->ip;
      }

      p2->iattrib[use_key] = k2.updateC(--C2);
      partflag     = 2;
      std::get<0>(ctd2->RR[id])++; 
      std::get<1>(ctd2->RR[id]) += NN;
      std::get<2>(ctd2->RR[id]) += kEe1[id] * NN;
    }

    delEeV = delE;

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
    double kEm = totE / dof;
    
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

int CollideIon::inelasticWeight(pHOT *tree, Particle* p1, Particle* p2, 
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
  
  if (ZWList[Z1] < ZWList[Z2]) {
    Particle *pT = p1;
    p1 = p2;
    p2 = pT;

    collTDPtr ctdT = ctd1;
    ctd1 = ctd2;
    ctd2 = ctdT;

    k1 = KeyConvert(p1->iattrib[use_key]);
    k2 = KeyConvert(p2->iattrib[use_key]);

    Z1 = k1.getKey().first;
    C1 = k1.getKey().second;
    Z2 = k2.getKey().first;
    C2 = k2.getKey().second;
  }
      
  // Find the trace ratio
  //
  double Wa = ZWList[Z1];
  double Wb = ZWList[Z2];
  double q  = Wb/Wa;

  // Number of atoms in each super particle
  //
  double N1 = p1->mass*UserTreeDSMC::Munit / amu * ZWList[Z1];
  double N2 = p2->mass*UserTreeDSMC::Munit / amu * ZWList[Z2];

  // Number of associated electrons for each particle
  //
  double ne1 = C1 - 1;
  double ne2 = C2 - 1;


  // The total mass in system units
  //
  double m1 = p1->mass * ZWList[Z1];
  double m2 = p2->mass * ZWList[Z2];

  double Mt = m1 + m2;
  if (Mt<=0.0) return ret;
  
  // Reduced mass in ballistic collision (system units)
  //
  double Mu = q * m1 * m2 / Mt;

  // Available center of mass energy in the ballistic collision
  // (system units)
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
    //          |
    //          v
    const bool DEBUG_F = false;
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
      std::get<1>(ctd1->ff[id]) += N1;
      std::get<2>(ctd1->ff[id]) += delE * N1;
    }

    if (interFlag == colexcite_1) {
      delE = IS.selectCEInteract(ch.IonList[Q1], CE1[id]);
      partflag      = 1;
      std::get<0>(ctd1->CE[id])++;
      std::get<1>(ctd1->CE[id]) += N1;
      std::get<2>(ctd1->CE[id]) += delE * N1;
    }

    if (interFlag == ionize_1) {
      delE          = IS.DIInterLoss(ch.IonList[Q1]);
      p1->iattrib[use_key] = k1.updateC(++C1);
      partflag      = 1;
      std::get<0>(ctd1->CI[id])++; 
      std::get<1>(ctd1->CI[id]) += N1;
      std::get<2>(ctd1->CI[id]) += delE * N1;
    }

    if (interFlag == recomb_1) {
      if (RECOMB_KE) {
	lQ rQ(Z1, C1-1);
	delE = ch.IonList[rQ]->ip + kEe1[id];
      }

      p1->iattrib[use_key] = k1.updateC(--C1);
      partflag      = 1;
      std::get<0>(ctd1->RR[id])++; 
      std::get<1>(ctd1->RR[id]) += N1;
      std::get<2>(ctd1->RR[id]) += kEe1[id] * N1;
    }
    
    //-------------------------
    // Particle 2 interactions
    //-------------------------

    if (interFlag == free_free_2) {
      delE          = IS.selectFFInteract(ch.IonList[Q2], id);
      partflag      = 2;
      std::get<0>(ctd2->ff[id])++;
      std::get<1>(ctd2->ff[id]) += N2;
      std::get<2>(ctd2->ff[id]) += delE * N2;
    }

    if (interFlag == colexcite_2) {
      delE         = IS.selectCEInteract(ch.IonList[Q2], CE2[id]);
      partflag     = 2;
      std::get<0>(ctd2->CE[id])++; 
      std::get<1>(ctd2->CE[id]) += N2;
      std::get<2>(ctd2->CE[id]) += delE * N2;
    }

    if (interFlag == ionize_2) {
      delE = IS.DIInterLoss(ch.IonList[Q2]);
      p2->iattrib[use_key] = k2.updateC(++C2);
      std::get<0>(ctd2->CI[id])++; 
      std::get<1>(ctd2->CI[id]) += N2;
      std::get<2>(ctd2->CI[id]) += delE * N2;
      partflag     = 2;
    }

    if (interFlag == recomb_2) {
      if (RECOMB_KE) {
	lQ rQ(Z1, C1-1);
	delE = ch.IonList[rQ]->ip + kEe2[id];
      }

      p2->iattrib[use_key] = k2.updateC(--C2);
      partflag     = 2;
      std::get<0>(ctd2->RR[id])++; 
      std::get<1>(ctd2->RR[id]) += N2;
      std::get<2>(ctd2->RR[id]) += kEe2[id] * N2;
    }

    delEeV = delE;

    // Convert to super particle
    //
    if (partflag == 1) delE *= N1;
    if (partflag == 2) delE *= N2;
    
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
  double remE=0.0, totE=kE, kEe=0.0, dof = 1.0;

  if (Z1 == Z2 and use_cons>=0) {
    totE = p1->dattrib[use_cons] + p2->dattrib[use_cons];
    p1->dattrib[use_cons] = p2->dattrib[use_cons] = 0.0;
  }

  if (use_Eint>=0) {
    totE += p1->dattrib[use_Eint] + p2->dattrib[use_Eint];
    p1->dattrib[use_Eint] = p2->dattrib[use_Eint] = 0.0;
  }


  // Electrons from Particle 2 have interacted with atom/ion in Particle 1
  //
  if (partflag==1) {
    remE = totE - dE;		// Energy floor
    kEe  = kEe2[id];		// Electron energy

				// Total degrees of freedom
    if (use_Eint>=0) dof = 1.0 + ne2;

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
    remE = totE - dE;		// Energy floor
    kEe  = kEe1[id];		// Electron energy

				// Total degrees of freedom
    if (use_Eint>=0) dof  = 1.0 + ne1;

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
    double kEm = totE / dof;
    
				// Get new relative velocity
    double ncr     = sqrt( 2.0*kEm/Mu );

    ret            = 0;		// No error

    if (partflag==1) {
      std::get<0>(ctd1->dv[id])++; 
      std::get<1>(ctd1->dv[id]) += N1;
      std::get<2>(ctd1->dv[id]) += 
	0.5*Mu*(vi - ncr)*(vi - ncr)/N1 * UserTreeDSMC::Eunit / eV;
    }
    
    if (partflag==2) {
      std::get<0>(ctd2->dv[id])++; 
      std::get<1>(ctd2->dv[id]) += N2;
      std::get<2>(ctd2->dv[id]) += 
	0.5*Mu*(vi - ncr)*(vi - ncr)/N2 * UserTreeDSMC::Eunit / eV;
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
				// All available energy will be lost
    lostSoFar[id] += totE;
    decolT[id]    += totE - delE;
    
    double ncr     = TolV;
    ret            = 1;		// Set error flag
    
				// Conservation of energy for internal
				// degrees of freedom
    totE          -= delE;

    if (partflag==1) {
      std::get<0>(ctd1->dv[id])++; 
      std::get<1>(ctd1->dv[id]) += N1;
      std::get<2>(ctd1->dv[id]) +=
	0.5*Mu*(vi - ncr)*(vi - ncr)/N1 * UserTreeDSMC::Eunit / eV;
    }
      
    if (partflag==2) {
      std::get<0>(ctd2->dv[id])++; 
      std::get<1>(ctd2->dv[id]) += N2;
      std::get<2>(ctd2->dv[id]) +=
	0.5*Mu*(vi - ncr)*(vi - ncr)/N2 * UserTreeDSMC::Eunit / eV;
    }
    
				// Reset internal energy excess
    if (use_exes>=0) {
      p1->dattrib[use_exes] = p1->mass*totE/Mt;
      p2->dattrib[use_exes] = p2->mass*totE/Mt;
    }
  }

  // Compute the new, scattered velocities
  //
  
  std::vector<double> vrel(3), vcm(3);

  for(unsigned k=0; k<3; k++)
    vcm[k] = (m1*p1->vel[k] + m2*p2->vel[k]) / Mt;
	    
  double cos_th = 1.0 - 2.0*(*unit)();       // Cosine and sine of
  double sin_th = sqrt(1.0 - cos_th*cos_th); // Collision angle theta
  double phi    = 2.0*M_PI*(*unit)();	     // Collision angle phi
  
  vrel[0] = vi*cos_th;	  // Compute post-collision relative
  vrel[1] = vi*sin_th*cos(phi); // velocity for an elastic 
  vrel[2] = vi*sin_th*sin(phi); // interaction

				// Center of mass lost energy
  double deltaKE = 0.0;
  for(unsigned k=0; k<3; k++) {
    double dv = -q*p1->vel[k] + q*(vcm[k] + p2->mass/Mt*vrel[k]);
    deltaKE += dv*dv;
  }
  deltaKE *= 0.5*ZWList[Z1]*p1->mass * q * (1.0 - q);
  
  if (Z1 == Z2)
    p1->dattrib[use_cons] = p2->dattrib[use_cons] = 0.5*deltaKE;
  else
    p1->dattrib[use_cons] = deltaKE;

  double vf = 0.0;
  if (totE > 0.0) {
    vf = sqrt( 2.0*totE/Mu );
  } else {
    if (Z1 == Z2)
      p1->dattrib[use_cons] = p2->dattrib[use_cons] = 0.5*totE;
    else
      p1->dattrib[use_cons] = totE;
  }

  double vfac = vf/vi;

  // Update post-collision velocities
  // 
  *cr = 0.0;
  for(unsigned k=0; k<3; k++ ) {
    p1->vel[k] = (1.0 - q)*p1->vel[k] + q*(vcm[k] + p2->mass/Mt*vrel[k]*vfac);
    p2->vel[k] = vcm[k] - p1->mass/Mt*vrel[k]*vfac;
    double dv = p1->vel[k] - p2->vel[k];
    *cr += dv*dv;
  }
  *cr = sqrt(*cr);
  
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

int CollideIon::inelasticTrace(pHOT *tree, Particle* p1, Particle* p2, 
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
  //          |
  //          v
  const bool DEBUG_F = false;
  //
  if (DEBUG_F) {
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
	  lQ rQ(Z, C-1);
	  double Xi = ch.IonList[rQ]->ip;
	  delE1 = (Xi + kEe2[id]) * N1;
	  delE2 = (Xi + kEe1[id]) * N2;
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

void CollideIon::finalize_cell(pHOT* tree, pCell* cell, double kedsp, int id)
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
	} else {
	  nOK = 1;		// Can't read use_cons value, fatal
	}

				// Print warning, not fatal
	if (use_cons<0) {
	  std::cout << "CollideIon: energy key position is not defined, "
		    << "you using trace-species weighting but not imposing energy conservation"
		    << std::endl;
	}

	if (nOK == 0) {
	  
	  int Z;
	  double W;
	  while (1) {
	    in.getline(line, nline);
	    if (in.good()) {
	      std::istringstream sz(line);
	      sz >> Z;
	      sz >> W;		// Add to the element list
	      if (!sz.bad()) {
		ZList.insert(Z);
		ZWList[Z] = W;
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

  } else {
    for (unsigned j=0; j<sz; j++) {
      MPI_Bcast(&z, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&v, 1, MPI_DOUBLE,         0, MPI_COMM_WORLD);
      ZWList[z] = v;
    }
  }

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
(pCell* c, sKeyDmap* Fn, double crm, double tau, int id,
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
(pCell* c, sKeyDmap* Fn, double crm, double tau, int id,
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
(pCell* c, sKeyDmap* Fn, double crm, double tau, int id,
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
    densM[i1] = c->Mass(i1) * ZWList[i1.first] / volc;
    //                        ^
    //                        |
    // Number density---------+
    //
  }
    
  if (0) {
    std::cout << std::setw(10) << "Species"
	      << std::setw(16) << "m density"
	      << std::setw(16) << "sp mass"
	      << std::setw(10) << "n count"
	      << std::setw(10) << "weight"
	      << std::endl
	      << std::setw(10) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(10) << "---------"
	      << std::setw(16) << "---------"
	      << std::endl;
    for (auto it : c->count) {
      std::ostringstream sout;
      sout << "(" << it.first.first << ", " << it.first.second << ")";
      std::cout << std::setw(10) << sout.str()
		<< std::setw(16) << densM[it.first]
		<< std::setw(16) << c->Mass(it.first)
		<< std::setw(10) << c->Count(it.first)
		<< std::setw(16) << ZWList[it.first.first]
		<< std::endl;
    }
  }

  double meanDens = 0.0;
  meanLambda      = 0.0;
  meanCollP       = 0.0;
    
  for (auto it1 : c->count) {

    if (it1.second == 0) continue;

    speciesKey i1 = it1.first;
    crossM [i1]   = 0.0;

    for (auto it2 : c->count) {

      if (it2.second == 0) continue;

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
	   << " (" << i1.first << ", " << i1.second << ")"
	   << " crossM = " << crossM[i1] 
	   << " densM = "  <<  densM[i1] 
	   << " Fn = "     <<  (*Fn)[i1] << endl;
      
      std::cout << std::setw(10) << "Species"
		<< std::setw(16) << "x-section"
		<< std::setw(16) << "m density"
		<< std::setw(16) << "sp mass"
		<< std::setw(10) << "n count"
		<< std::endl
		<< std::setw(10) << "---------"
		<< std::setw(16) << "---------"
		<< std::setw(16) << "---------"
		<< std::setw(16) << "---------"
		<< std::setw(10) << "---------"
		<< std::endl;
      for (auto it : csections[id][i1]) {
	std::ostringstream sout;
	sout << "(" << it.first.first << ", " << it.first.second << ")";
	cout << std::setw(10) << sout.str()
	     << std::setw(16) << it.second 
	     << std::setw(16) << densM[it.first]
	     << std::setw(16) << c->Mass(it.first)
	     << std::setw(10) << c->Count(it.first)
	     << std::endl;
      }
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
      double Prob = densM[i1] * (*Fn)[i2] * csections[id][i1][i2] * crm * tau;
      
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

    double mass = 0.0;
    tempM = 0.0;

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
    }

    // Send values to root
    //
    double val1, val2;
    
    for (int i=1; i<numprocs; i++) {
      if (i == myid) {
				// Mass
	MPI_Send(&mass,  1, MPI_DOUBLE, 0, 331, MPI_COMM_WORLD);
				// Temp
	MPI_Send(&tempM, 1, MPI_DOUBLE, 0, 332, MPI_COMM_WORLD);
      }

				// Root receives from Node i
      if (0 == myid) {
	MPI_Recv(&val1, 1, MPI_DOUBLE, i, 331, MPI_COMM_WORLD, 
		 MPI_STATUS_IGNORE);
	MPI_Recv(&val2, 1, MPI_DOUBLE, i, 332, MPI_COMM_WORLD, 
		 MPI_STATUS_IGNORE);
	mass  += val1;
	tempM += val2;
      }
    }

    if (mass>0.0) tempM /= mass;
    
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

  if (aType==Direct) {
    double numbC = 0.0, massC = 0.0;
    for (auto it : cell->count) {
      speciesKey i = it.first;
      double M = cell->Mass(i);
      numbC += M / atomic_weights[i.first];
      massC += M;
    }

    mol_weight = massC/numbC;
  }

  if (aType==Weight) {
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
