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

// Proton mass (g)
const double mp      = 1.67262158e-24;

// Electron mass (g)
const double me      =  9.10938291e-28;

// Speed of light (cm/s)
const double c       =  2.99792458e+10;

// electron volt in (cgs)
const double eV      =  1.60217653e-12;

// Boltzmann constant (cgs)
const double boltz   = 1.3810e-16;

// Boltzmann constant (eV)
const double boltzEv = 8.6173324e-5;

// Planck's constant (cgs)
//const double planck  = 6.6262e-27;

// Electron charge (cgs)
const double esu     = 4.80320427e-10;

// Atomic mass unit in grams
const double amu     = 1.660539e-24;

// Number of Ion elements
const int N_Z = 2;		

// Array of the different Z's for the elements
const int ZList[N_Z] = {1, 2}; 

double   CollideIon::Nmin    = 1.0e-08;
double   CollideIon::Nmax    = 1.0e+25;
double   CollideIon::Tmin    = 1.0e+03;
double   CollideIon::Tmax    = 1.0e+08;
double   CollideIon::TolV    = 1.0e-03;
unsigned CollideIon::Nnum    = 400;
unsigned CollideIon::Tnum    = 200;
string   CollideIon::cache   = ".HeatCool";
bool     CollideIon::frost_warning = false;

bool NO_COOL = false;

// Minimum energy for Rutherford scattering of ions
const double FloorEv = 0.05;

CollideIon::CollideIon(ExternalForce *force, double hD, double sD, int Nth) : 
  Collide(force, hD, sD, Nth)
{
  NUM = 0;
  csections = std::vector<sKey2Dmap> (nthrds);
  
  for(int i = 0; i < N_Z; i++) {
    for (int j = 1; j <= ZList[i] + 1; j++) {
      IonList[ZList[i]][j] = Ion(ZList[i], j, ch);
      IonList[ZList[i]][j].freeFreeDifferential(ch);
    }
  }

  for (int i = 0; i < N_Z; i++) {
    Ni.push_back(1);		// dummy value for now
  }

  // Random variable generators
  gen  = new ACG(11+myid);
  unit = new Uniform(0.0, 1.0, gen);
  
  // Energy diagnostics
  totalSoFar = 0.0;
  massSoFar  = 0.0;
  lostSoFar  = vector<double>(nthrds, 0.0);
  
  /** 
      Initialize the atomic_weights map
      Hardcode the atomic weight map for use in collFrac
  */
  atomic_weights[1]  = 1.0079;
  atomic_weights[2]  = 4.0026;
  atomic_weights[3]  = 6.941;
  atomic_weights[4]  = 9.0122;
  atomic_weights[5]  = 10.811;
  atomic_weights[6]  = 12.011;
  atomic_weights[7]  = 14.007;
  atomic_weights[8]  = 15.999;
  atomic_weights[9]  = 18.998;
  atomic_weights[10] = 20.180;
  atomic_weights[11] = 22.990;
  atomic_weights[12] = 24.305;
  
  collD = boost::shared_ptr<collDiag>(new collDiag());

  if (myid==0 && NO_COOL) std::cout << "No cooling is set to TRUE" 
				    << std::endl;
}

CollideIon::~CollideIon()
{
  /** I don't believe this is needed.  Destroying the vector will call
      the destructor for Ion. */
  IonList.erase(IonList.begin(), IonList.end());
}


/**
   Precompute all the necessary cross sections
 */
void CollideIon::initialize_cell(pHOT* tree, pCell* cell, double rvmax, int id)
{
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
  double EeV = Eerg / eV;

  // Upscaling factor for scattering cross section
  //
  double sUp  = diamfac*diamfac;

  // Electron velocity equipartition factor
  double eVel = sqrt(mp/me);

  typedef std::map<speciesKey, unsigned> Count;

  for (Count::iterator it1 = cell->count.begin(); it1 != cell->count.end(); it1++)  {

    speciesKey i1 = it1->first;
    double Cross1 = geometric(i1.first);
    
    for (Count::iterator it2 = cell->count.begin(); it2 != cell->count.end(); it2++) {
      
      speciesKey i2 = it2->first;
      double Cross2 = geometric(i2.first);

      double mu = atomic_weights[i1.first] * atomic_weights[i1.first] / 
	(atomic_weights[i1.first] + atomic_weights[i2.first]);

      if (i2.second>1) {
	double ne2 = i2.second - 1;
	if (i1.second==1)
	  Cross1 = elastic(i1.first, EeV * mu) * eVel * ne2;
	else {
	  double b = 0.5*esu*esu*(i1.second - 1) /
	    std::max<double>(Eerg, FloorEv*eV) * 1.0e7; // nm
	  Cross1 = M_PI*b*b * eVel * ne2;
	}
      }

      if (i1.second>1) {
	double ne1 = i1.second - 1;
	if (i2.second==1)
	  Cross2 = elastic(i2.first, EeV * mu) * eVel * ne1;
	else {
	  double b = 0.5*esu*esu*(i2.second - 1) /
	    std::max<double>(Eerg, FloorEv*eV) * 1.0e7; // nm
	  Cross2 = M_PI*b*b * eVel * ne1;
	}
      }

      
      csections[id][i1][i2] = (Cross1 + Cross2) * sUp * 1e-14 / 
	(UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
    }
  }
}


sKey2Dmap& CollideIon::totalScatteringCrossSections(double crm, pCell *c, int id)
{
  typedef std::map<speciesKey, unsigned> Count;
  Count::iterator it1, it2;
  
  double vel = crm * UserTreeDSMC::Vunit;
  double Eerg = 0.5*vel*vel*amu;
  double EeV  = Eerg / eV;

  // Upscaling factor for scattering cross section
  //
  double sUp  = diamfac*diamfac;

  for (it1 = c->count.begin(); it1 != c->count.end(); it1++)  {

    speciesKey i1 = it1->first;
    double Cross1 = geometric(i1.first);
    
    for (it2 = c->count.begin(); it2 != c->count.end(); it2++)  
      {
	speciesKey i2 = it2->first;

	double Cross2 = geometric(i2.first);

	double mu = atomic_weights[i1.first] * atomic_weights[i1.first] / 
	  (atomic_weights[i1.first] + atomic_weights[i2.first]);

	double eVel1 = sqrt(amu*atomic_weights[i1.first]/me);
	double eVel2 = sqrt(amu*atomic_weights[i2.first]/me);

	// Electrons in second particle
	//
	if (i2.second>1) {
	  double ne2 = i2.second - 1;
	  if (i1.second==1)	// Neutral atom-electron scattering
	    Cross1 = elastic(i1.first, EeV * mu) * eVel2*ne2;
	  else {		// Rutherford scattering
	    double b = 0.5*esu*esu*(i1.second - 1) /
	      std::max<double>(Eerg*mu, FloorEv*eV) * 1.0e7; // nm
	    Cross1 = M_PI*b*b * eVel2*ne2;
	  }
	}

	// Electrons in first particle
	//
	if (i1.second>1) {
	  double ne1 = i1.second - 1;
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
    
  return csections[id];
}

double CollideIon::crossSection(pHOT *tree, Particle* p1, Particle* p2, 
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
  kEi = 0.5 * mu * vel*vel;

  // Electron velocity equipartition factors
  //
  double eVel1 = sqrt(m1/me);
  double eVel2 = sqrt(m2/me);

  // Internal energy per particle
  //
  Ein1 = Ein2 = 0.0;
  if (use_Eint>=0) {
    Ein1 = p1->dattrib[use_Eint] * UserTreeDSMC::Eunit/N1;
    Ein2 = p2->dattrib[use_Eint] * UserTreeDSMC::Eunit/N2;
  }

  // Compute the total available energy and divide among degrees of freedom
  // Convert ergs to eV
  //

  kEe1 = (kEi + Ein1)/(1.0 + ne1) / eV;
  kEe2 = (kEi + Ein2)/(1.0 + ne2) / eV;
  
  // Get temperatures from cells
  //
  double E_therm_1 = kEe1;
  double E_therm_2 = kEe2;

  // Or use the effective temperature from the COM KE
  //
  if (0 && use_temp>=0) {
    if (p1->dattrib[use_temp] > 0) E_therm_1 = boltzEv * p1->dattrib[use_temp];
    if (p2->dattrib[use_temp] > 0) E_therm_2 = boltzEv * p2->dattrib[use_temp];
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

  double cross12 = geometric(Z1);
  double cross21 = geometric(Z2);
	
  // Electrons in second particle
  //
  if (ne2 > 0) {
    if (C1==1)		// Neutral atom-electron scattering
      cross12 = elastic(Z1, kEe2) * eVel2*ne2;
    else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C1-1) /
	std::max<double>(kEe2*eV, FloorEv*eV) * 1.0e7; // nm
      cross12 = M_PI*b*b * eVel2*ne2;
    }
  }
    
  // Electrons in first particle
  //
  if (ne1 > 0) {
    if (C2==1)		// Neutral atom-electron scattering
      cross21 = elastic(Z2, kEe1) * eVel1*ne1;
    else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C2-1) /
	std::max<double>(kEe1*eV, FloorEv*eV) * 1.0e7; // nm
      cross21 = M_PI*b*b * eVel1*ne1;
    }
  }

  dCrossMap[id].push_back((cross12 + cross21)*diamfac*diamfac);
  dInterMap[id].push_back(0);


  //--------------------------------------------------
  // Particle 1 interacts with Particle 2
  //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
  if (C1 > 1 and ne2 > 0) {
    double ff1 = IonList[Z1][C1].freeFreeCross(ch, kEe2);
    double crs = eVel2*ne2 * ff1;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(1);

    sum12 += crs;
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  if (ne2 > 0 and C1 <= Z1) {

    if (E_therm_2 == 0.0) {
      std::cout << "E_therm_2 == 0.0, kEe=" << kEe2 << ", cr=" 
		<< cr << std::endl;
    }
    
    CE1[id] = IonList[Z1][C1].collExciteCross(ch, kEe2, E_therm_2);
    double crs = eVel2*ne2 * CE1[id].back().first;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(2);
    sum12 += crs;
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  if (C1 < (Z1 + 1) and ne2 > 0) {

    double DI1 = IonList[Z1][C1].directIonCross(ch, kEe2);
    double crs = eVel2*ne2 * DI1;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(3);

    sum12 += crs;
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  if (C1 > 1 and ne2 > 0) {

    std::vector<double> RE1 = IonList[Z1][C1].radRecombCross(ch, kEe2);
    double crs = eVel2*ne2 * RE1.back();

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(4);

    sum12 += crs;
  }

  
  //--------------------------------------------------
  // Particle 2 interacts with Particle 1
  //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
  if (C2 > 1 and ne1 > 0) {
    double ff2 = IonList[Z2][C2].freeFreeCross(ch, kEe1);
    double crs = eVel1*ne1 * ff2;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(6);

    sum21 += crs;
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  if (ne1 > 0 and C2 <= Z2) {

    if (E_therm_1 == 0.0) {
      std::cout << "E_therm_1 == 0.0, kEe1=" << kEe1 << ", cr=" << cr 
		<< std::endl;
    }

    CE2[id] = IonList[Z2][C2].collExciteCross(ch, kEe1, E_therm_1);
    double crs = eVel1*ne1 * CE2[id].back().first;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(7);

    sum21 += crs;
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  if (C2 < (Z2 + 1) and ne1 > 0) {
    double DI2 = IonList[Z2][C2].directIonCross(ch, kEe1);
    double crs = ne1 * DI2;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(8);

    sum21 += crs;
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  if (C2 > 1 and ne1 > 0) {
    std::vector<double> RE2 = IonList[Z2][C2].radRecombCross(ch, kEe1);
    double crs = eVel1*ne1*RE2.back();

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(9);
    
    sum21 += crs;
  } 
  
				//-------------------------------
				// *** Convert to system units
				//-------------------------------
  return (cross12 + cross21 + sum12 + sum21) * 1e-14 / 
    (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
}


int CollideIon::inelastic(pHOT *tree, Particle* p1, Particle* p2, 
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
    // Particle 1 interactions
    //-------------------------

    if (interFlag == 1) {
      delE          = IS.selectFFInteract(IonList[Z1][C1], kEe2);
      partflag      = 1;
      ctd1->ff[id].first++; 
      ctd1->ff[id].second  += delE;
    }

    if (interFlag == 2) {
      delE = IS.selectCEInteract(IonList[Z1][C1], CE1[id]);
      partflag      = 1;
      ctd1->CE[id].first++; 
      ctd1->CE[id].second  += delE;
    }

    if (interFlag == 3) {
      delE          = IS.DIInterLoss(ch, IonList[Z1][C1]);
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
    if (interFlag == 4) {
      delE          = kEe2;
      p1->iattrib[use_key] = k1.updateC(--C1);
      assert(C1 > 0);
      partflag      = 1;
      ctd1->RR[id].first++; 
      ctd1->RR[id].second  += delE;
    }
    
    //-------------------------
    // Particle 2 interactions
    //-------------------------

    if (interFlag == 6) {
      delE          = IS.selectFFInteract(IonList[Z2][C2], kEe1);
      partflag      = 2;
      ctd2->ff[id].first++;
      ctd2->ff[id].second += delE;
    }

    if (interFlag == 7) {
      delE         = IS.selectCEInteract(IonList[Z2][C2], CE2[id]);
      partflag     = 2;
      ctd2->CE[id].first++; 
      ctd2->CE[id].second += delE;
    }

    if (interFlag == 8) {
      delE = IS.DIInterLoss(ch, IonList[Z2][C2]);
      p2->iattrib[use_key] = k2.updateC(++C2);
      ctd2->CI[id].first++; 
      ctd2->CI[id].second += delE;
      partflag     = 2;
    }

    if (interFlag == 9) {
      delE         = kEe1;	// See comment above for interFlag==4
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
    delE = delE*1.602177e-12;
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
  double remE, totE, kEe, dof;

  // Electrons from Particle 2 have interacted with atom/ion in Particle 1
  //
  if (partflag==1) {
    totE = kE;			// KE + internal
    if (use_Eint>=0) totE += p2->dattrib[use_Eint];
    remE = totE - dE;		// Energy floor
    kEe  = kEe2;		// Electron energy
    dof  = 1.0 + ne2;		// Total degrees of freedom

				// Energy diagnostics
    ctd1->eV_av[id] += kEe2;
    if (std::isnan(ctd2->eV_av[id])) {
      std::cout << "eV_N=" << ctd1->eV_N[id] << std::endl;
    }
    ctd1->eV_N[id]++;
    ctd1->eV_min[id] = std::min(ctd1->eV_min[id], kEe2);
    ctd1->eV_max[id] = std::max(ctd2->eV_max[id], kEe2);
    
    if (kEe2 > 10.2) { ctd1->eV_10[id]++;}
  }

  // Electrons from Particle 1 interacted with atom/ion in Particle 2
  //
  if (partflag==2) {
    totE = kE;			// KE + internal
    if (use_Eint>=0) totE += p1->dattrib[use_Eint];
    remE = totE - dE;		// Energy floor
    kEe  = kEe1;		// Electron energy
    dof  = 1.0 + ne1;		// Total degrees of freedom

				// Energy diagnostics
    ctd2->eV_av[id] += kEe1;
    if (std::isnan(ctd2->eV_av[id])) {
      std::cout << "eV_N=" << ctd2->eV_N[id] << std::endl;
    }
    ctd2->eV_N[id]++;
    ctd2->eV_min[id] = std::min(ctd2->eV_min[id], kEe1);
    ctd2->eV_max[id] = std::max(ctd2->eV_max[id], kEe1);
    
    if (kEe1 > 10.2) { ctd2->eV_10[id]++; }
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
	0.5*Mu*(vi - (*cr))*(vi - (*cr)) * UserTreeDSMC::Eunit / eV;
    }
    
    if (partflag==2) {
      ctd2->dv[id].first++; 
      ctd2->dv[id].second += 
	0.5*Mu*(vi - (*cr))*(vi - (*cr)) * UserTreeDSMC::Eunit / eV;
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
	  0.5*Mu*(vi - (*cr))*(vi - (*cr)) * UserTreeDSMC::Eunit / eV;
      }
      
      if (partflag==2) {
	ctd2->dv[id].first++; 
	ctd2->dv[id].second +=
	  0.5*Mu*(vi - (*cr))*(vi - (*cr)) * UserTreeDSMC::Eunit / eV;
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
	  0.5*Mu*(vi - (*cr))*(vi - (*cr)) * UserTreeDSMC::Eunit / eV;
      }
      
      if (partflag==2) {
	ctd2->dv[id].first++; 
	ctd2->dv[id].second +=
	  0.5*Mu*(vi - (*cr))*(vi - (*cr)) * UserTreeDSMC::Eunit / eV;
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
  
  // thread_timing_beg(id);
  
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
  
  //thread_timing_end(id);
  
  return (NULL);
}

void CollideIon::finalize_cell(pHOT* tree, pCell* cell, double kedsp, int id)
{
  
}

collDiag::collDiag()
{
  // Initialize the map
  //
  for (int n=0; n<N_Z; n++) {
    unsigned short Z = ZList[n];
    for (unsigned short C=1; C<Z+2; C++) {
      speciesKey k(Z, C);
      (*this)[k] = collTDPtr(new CollisionTypeDiag());
    }
  }

  initialize();
}

void collDiag::gather()
{
  for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
    collTDPtr ctd = it->second;
    ctd->sumUp();

    unsigned u;
    double z;
    if (myid==0) {
      MPI_Reduce(&(u=ctd->ff_s.first),  &ctd->ff_s.first,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 1
      MPI_Reduce(&(z=ctd->ff_s.second), &ctd->ff_s.second, 1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 2
      MPI_Reduce(&(u=ctd->CE_s.first),  &ctd->CE_s.first,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 3
      MPI_Reduce(&(z=ctd->CE_s.second), &ctd->CE_s.second, 1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 4
      MPI_Reduce(&(u=ctd->CI_s.first),  &ctd->CI_s.first,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 5
      MPI_Reduce(&(z=ctd->CI_s.second), &ctd->CI_s.second, 1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 6
      MPI_Reduce(&(u=ctd->RR_s.first),  &ctd->RR_s.first,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 7
      MPI_Reduce(&(z=ctd->RR_s.second), &ctd->RR_s.second, 1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 8
      MPI_Reduce(&(u=ctd->dv_s.first),  &ctd->dv_s.first,  1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 9
      MPI_Reduce(&(z=ctd->dv_s.second), &ctd->dv_s.second, 1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 10
      MPI_Reduce(&(z=ctd->eV_av_s),     &ctd->eV_av_s,     1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 11
      MPI_Reduce(&(u=ctd->eV_N_s),      &ctd->eV_N_s,      1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 12
      MPI_Reduce(&(z=ctd->eV_min_s),    &ctd->eV_min_s,    1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 13
      MPI_Reduce(&(z=ctd->eV_max_s),    &ctd->eV_max_s,    1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 14
      MPI_Reduce(&(u=ctd->eV_10_s),     &ctd->eV_10_s,     1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 15
    } else {
      MPI_Reduce(&ctd->ff_s.first,      &u,                1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 1 
      MPI_Reduce(&ctd->ff_s.second, 	&z, 		   1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 2 
      MPI_Reduce(&ctd->CE_s.first,  	&u, 		   1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 3 
      MPI_Reduce(&ctd->CE_s.second, 	&z, 		   1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 4 
      MPI_Reduce(&ctd->CI_s.first,  	&u, 		   1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 5 
      MPI_Reduce(&ctd->CI_s.second, 	&z, 		   1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 6 
      MPI_Reduce(&ctd->RR_s.first,  	&u, 		   1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 7 
      MPI_Reduce(&ctd->RR_s.second, 	&z, 		   1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 8 
      MPI_Reduce(&ctd->dv_s.first,    	&u, 		   1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 9 
      MPI_Reduce(&ctd->dv_s.second,   	&z, 		   1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 10
      MPI_Reduce(&ctd->eV_av_s,       	&z, 		   1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 11
      MPI_Reduce(&ctd->eV_N_s,        	&u, 		   1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 12
      MPI_Reduce(&ctd->eV_min_s,      	&z, 		   1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 13
      MPI_Reduce(&ctd->eV_max_s,      	&z, 		   1, MPI_DOUBLE,   MPI_SUM, 0, MPI_COMM_WORLD); // 14
      MPI_Reduce(&ctd->eV_10_s,       	&u, 		   1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD); // 15
    }
  }
}

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
	    << "#"                                         << std::endl;
	
	out << "#" << std::setw(11) << "Species==>" << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  ostringstream sout, sout2;
	  sout  << "(" << it->first.first << ", " << it->first.second << ")";
	  size_t w = 9*12, l = sout.str().size();
	  sout2 << std::setw((w-l)/2) << ' ' << sout.str();
	  out   << std::setw(w) << sout2.str() << " | ";
	}
	out << std::endl;

	out << "#" << std::setw(11) << "Time" << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  out << std::setw(12) << "N(ff) "
	      << std::setw(12) << "E(ff) "
	      << std::setw(12) << "N(ce) "
	      << std::setw(12) << "E(ce) "
	      << std::setw(12) << "N(ci) "
	      << std::setw(12) << "E(ci) "
	      << std::setw(12) << "N(rr) "
	      << std::setw(12) << "E(rr) "
	      << std::setw(12) << "d(KE) "
	      << " | ";
	}
	out << std::endl;
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

	out << "#" << std::setw(11) << "Species==>" << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  ostringstream sout, sout2;
	  sout  << "(" << it->first.first << ", " << it->first.second << ")";
	  size_t w = 5*12, l = sout.str().size();
	  sout2 << std::setw((w-l)/2) << ' ' << sout.str();
	  out   << std::setw(w) << sout2.str() << " | ";
	}
	out << std::endl;
	
	out << "#" << std::setw(11) << "Time" << " | ";
	for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	  out << std::setw(12) << "avg "
	      << std::setw(12) << "num "
	      << std::setw(12) << "min "
	      << std::setw(12) << "max "
	      << std::setw(12) << "over10 "
	      << " | ";
	}
	out << std::endl;
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
    if (out) {
      out << std::setw(12) << tnow << " | ";
      for (sKeyCollTD::iterator it=this->begin(); it!=this->end(); it++) {
	collTDPtr ctd = it->second;
	out << std::setw(12) << ctd->ff_s.first
	    << std::setw(12) << ctd->ff_s.second
	    << std::setw(12) << ctd->CE_s.first
	    << std::setw(12) << ctd->CE_s.second
	    << std::setw(12) << ctd->CI_s.first
	    << std::setw(12) << ctd->CI_s.second
	    << std::setw(12) << ctd->RR_s.first
	    << std::setw(12) << ctd->RR_s.second;
	if (ctd->dv_s.first>0.0)
	  out << std::setw(12) << ctd->dv_s.second/ctd->dv_s.first << " | ";
	else
	  out << std::setw(12) << 0.0 << " | ";
      }
      out << std::endl;
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

