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
  
  // Make sure that internal energy is assigned
  //
  if (use_Eint<0) {
    if (myid==0) {
      std::cout << "*** Internal energy variable not assigned; "
		<< "will translational energy only, but this "
		<< " is probably not what you intend ***"
		<< std::endl;
    }
  }

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
  
  ff_d   = std::make_pair(0., 0.);
  CE_d   = std::make_pair(0., 0.);
  CI_d   = std::make_pair(0., 0.);
  RR_d   = std::make_pair(0., 0.);
  dv     = std::make_pair(0., 0.);

  eV_av  = 0;
  eV_N   = 0;
  eV_min = 999999;
  eV_max = 0;
  eV_10  = 0;
  
  printCollInitialize();

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

  // Electron velocity equipartition factor
  //
  double eVel = sqrt(mu/me);

  // Internal energy per particle
  //
  double dof = 1.0;
  Ein = 0.0;
  if (use_Eint>=0) {
    Ein += p1->dattrib[use_Eint] * UserTreeDSMC::Eunit/N1;
    Ein += p2->dattrib[use_Eint] * UserTreeDSMC::Eunit/N2;
    dof += ne1 + ne2;
  }

  // Compute the total available energy and divide among degrees of freedom
  // Convert ergs to eV
  //

  kEe = (kEi + Ein)/dof / eV;

  // Get temperatures from cells
  //
  double E_therm_1 = boltzEv * p1->dattrib[use_temp];
  double E_therm_2 = boltzEv * p2->dattrib[use_temp];
  
  // Or use the effective temperature from the COM KE
  //
  if (p1->dattrib[use_temp] == 0) {
    E_therm_1 = kEe;
  }

  if (p2->dattrib[use_temp] == 0) {
    E_therm_2 = kEe;
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
      cross12 = elastic(Z1, kEe) * eVel*ne2;
    else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C1-1) /
	std::max<double>(kEe*eV, FloorEv*eV) * 1.0e7; // nm
      cross12 = M_PI*b*b * eVel*ne2;
    }
  }
    
  // Electrons in first particle
  //
  if (ne1 > 0) {
    if (C2==1)		// Neutral atom-electron scattering
      cross21 = elastic(Z2, kEe) * eVel*ne1;
    else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C2-1) /
	std::max<double>(kEe*eV, FloorEv*eV) * 1.0e7; // nm
      cross21 = M_PI*b*b * eVel*ne1;
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
    double ff1 = IonList[Z1][C1].freeFreeCross(ch, kEe);
    double crs = eVel*ne2 * ff1;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(1);

    sum12 += crs;
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  if (ne2 > 0 and C1 <= Z1) {

    if (E_therm_1 == 0.0) {
      std::cout << "E_therm_1 == 0.0, kEe=" << kEe << ", cr=" << cr << std::endl;
    }
    
    CE1[id] = IonList[Z1][C1].collExciteCross(ch, kEe, E_therm_1);
    double crs = eVel*ne2 * CE1[id].back().first;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(2);
    sum12 += crs;
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  double DI1;
  if (C1 < (Z1 + 1) and ne2 > 0) {

    DI1 = IonList[Z1][C1].directIonCross(ch, kEe);
    double crs = eVel*ne2 * DI1;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(3);

    sum12 += crs;
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  std::vector<double> RE1;
  if (C1 > 1 and ne2 > 0) {

    RE1 = IonList[Z1][C1].radRecombCross(ch, kEe);
    double crs = eVel*ne2 * RE1.back();

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
    double ff2 = IonList[Z2][C2].freeFreeCross(ch, kEe);
    double crs = eVel*ne1 * ff2;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(6);

    sum21 += crs;
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  if (ne1 > 0 and C2 <= Z2) {

    if (E_therm_2 == 0.0) {
      std::cout << "E_therm_2 == 0.0, kEe=" << kEe << ", cr=" << cr 
		<< std::endl;
    }

    CE2[id] = IonList[Z2][C2].collExciteCross(ch, kEe, E_therm_2);
    double crs = eVel*ne1 * CE2[id].back().first;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(7);

    sum21 += crs;
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  double DI2;
  if (C2 < (Z2 + 1) and ne1 > 0) {
    DI2 = IonList[Z2][C2].directIonCross(ch, kEe);
    double crs = ne1 * DI2;

    dCrossMap[id].push_back(crs);
    dInterMap[id].push_back(8);

    sum21 += crs;
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  std::vector<double> RE2;
  if (C2 > 1 and ne1 > 0) {
    RE2 = IonList[Z2][C2].radRecombCross(ch, kEe);
    double crs = eVel*ne1*RE2.back();

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

  unsigned short Z1 = k1.getKey().first, C1 = k1.getKey().second;
  unsigned short Z2 = k2.getKey().first, C2 = k2.getKey().second;
  
  // Number of atoms in each super particle
  //
  double N1 = (p1->mass*UserTreeDSMC::Munit)/(atomic_weights[Z1]*amu);
  double N2 = (p2->mass*UserTreeDSMC::Munit)/(atomic_weights[Z2]*amu);
  double NN = std::min<double>(N1, N2);
  
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
  double dof  = 1.0;
				// KE + internal
  double totE = kE;
  if (use_Eint>=0) {
    dof  += ne1 + ne2;
    totE += p1->dattrib[use_Eint] + p2->dattrib[use_Eint];
  }

  //
  // Remainder energy after removing floor (system units)
  double remE = totE - dE;
  double delE = 0.0, delEeV = 0.0;

  // Now that the interactions have been calculated, create the
  // normalized cross section list to pick the interaction
  //
  std::vector<double> TotalCross;
  double tCross;
  tCross = 0;
  int si = 0;
  for (size_t i = 0; i < dCrossMap[id].size(); i++) {
    tCross += dCrossMap[id][i];
    TotalCross.push_back(tCross);
    si++;
  }

  assert (TotalCross.size() == dCrossMap[id].size());

  if (tCross != 0) {
    std::vector<double> CDF;
    std::vector <double>::iterator 
      max = max_element(TotalCross.begin(), TotalCross.end());

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

    // Sanity check: did not set interaction type??
    //
    if (interFlag<0) {
      std::cout << "interFlag NOT set, myid=" << myid 
		<< ", index=" << index << std::endl;
      index = 0;
      interFlag = dInterMap[id][0];
    }

    int partflag = 0;		// Will be 1 or 2, dependending on
				// which ion or neutral is selected
				// for interaction

    //-------------------------
    // Particle 1 interactions
    //-------------------------

    if (interFlag == 1) {
      delE          = IS.selectFFInteract(IonList[Z1][C1], kEe);
      partflag      = 1;
      ff_d.first++; 
      ff_d.second  += delE;
    }

    if (interFlag == 2) {
      delE = IS.selectCEInteract(IonList[Z1][C1], CE1[id]);
      partflag      = 1;
      CE_d.first++; 
      CE_d.second  += delE;
    }

    if (interFlag == 3) {
      delE          = IS.DIInterLoss(ch, IonList[Z1][C1]);
      C1++;
      assert(C1 <= (Z1 + 1));
      partflag      = 1;
      CI_d.first++; 
      CI_d.second  += delE;
    }

    // KE carried by electron is subtracted from the thermal reservoir
    // but not the radiation of "binding".  This radiation decreases
    // the total energy of the gas but not the thermal component.
    //
    if (interFlag == 4) {
      delE          = kEe;
      C1--;
      assert(C1 > 0);
      partflag      = 1;
      RR_d.first++; 
      RR_d.second  += delE;
    }
    
    //-------------------------
    // Particle 2 interactions
    //-------------------------

    if (interFlag == 6) {
      delE          = IS.selectFFInteract(IonList[Z2][C2], kEe);
      partflag      = 2;
      ff_d.first++;
      ff_d.second += delE;
    }

    if (interFlag == 7) {
      delE         = IS.selectCEInteract(IonList[Z2][C2], CE2[id]);
      partflag     = 2;
      CE_d.first++; 
      CE_d.second += delE;
    }

    if (interFlag == 8) {
      delE = IS.DIInterLoss(ch, IonList[Z2][C2]);
      C2++;
      CI_d.first++; 
      CI_d.second += delE;
      partflag     = 2;
    }

    if (interFlag == 9) {
      delE         = kEe;	// See comment above for interFlag==4
      C2--;
      assert(C2 > 0);
      partflag     = 2;
      RR_d.first++; 
      RR_d.second += delE;
    }

    delEeV = delE;

    // Convert to super particle
    //
    if (partflag) delE *= NN;
    
    // Energy summary diagnostics
    //
    eV_av += kEe;
    if (std::isnan(eV_av)) {
      std::cout << "eV_N=" << eV_N << std::endl;
    }
    eV_N++;
    eV_min = std::min(eV_min, kEe);
    eV_max = std::max(eV_max, kEe);
    
    if (kEe > 10.2) { eV_10++;}
  

    // Convert back to cgs
    //
    delE = delE*1.602177e-12;
  }
  
  assert(delE >= 0.0);

  // Artifically prevent cooling
  //
  if (NO_COOL) delE = 0.0;

  // Convert energy loss to system units
  //
  delE = delE/UserTreeDSMC::Eunit;
  
  
  if (frost_warning && delE*0.999 > totE) {
    std::cout << "delE > KE!! (" << delE << " > " << totE
	      << "), Interaction type = " << interFlag 
	      << " kEe  = "  << kEe
	      << " delE = " << delEeV 
	      << std::endl;
  }
  
  // Distribute electron energy among particles in pure scattering event
  //
  if (delE<=0.0) {
    if (use_Eint>=0 && ne1 + ne2>0) {
      double kEm  = (p1->dattrib[use_Eint] + p2->dattrib[use_Eint])/(ne1 + ne2);
      p1->dattrib[use_Eint] = ne1 * kEm;
      p2->dattrib[use_Eint] = ne2 * kEm;
    }
    return ret;
  }
  
  // Cooling rate diagnostic
  //
  if (TSDIAG) {
    if (delE>0.0) {
      int indx = (int)floor(log(remE/delE)/(log(2.0)*TSPOW) + 5);
      if (indx<0 ) indx = 0;
      if (indx>10) indx = 10;
      
      EoverT[id][indx] += Mt;
    }
  }
  
  if (use_exes>=0) {
    // (-/+) value means under/overcooled: positive/negative increment
    // to delE NB: delE may be < 0 if too much energy was radiated
    // previously . . .
    //
    delE -= p1->dattrib[use_exes] + p2->dattrib[use_exes];
  }
  
  // Sufficient energy available
  //
  if (remE > delE) {
    double vi      = (*cr);
    lostSoFar[id] += delE;
    decelT[id]    += delE;

    totE          -= delE;	// Remove the energy from the total
				// available

    double kEm     = totE/dof;	// Energy per particle

				// Get new relative velocity
    (*cr)          = sqrt( 2.0*kEm/Mu );

    dv.first++; 
    dv.second     += 0.5*(vi - (*cr))*(vi - (*cr));
    ret            = 0;		// No error

				// Distribute electron energy to particles
    if (use_Eint>=0) {
      p1->dattrib[use_Eint] = ne1 * kEm;
      p2->dattrib[use_Eint] = ne2 * kEm;
    }
				// Zero out internal energy excess
    if (use_exes>=0)		// since excess is now used up
      p1->dattrib[use_exes] = p2->dattrib[use_exes] = 0.0;

  } else {			// Inconsistent: too much energy lost!

    if (frost_warning && remE*0.999 > delE) {
      std::cout << "************ remE < delE!!!!" 
		<< " RemE = " << remE << " delE = " << delE 
		<< " interaction: " << interFlag    << std::endl;
    }
    
    lostSoFar[id] += remE;
    decolT[id]    += remE - delE;

    (*cr)         *= TolV;
    ret            = 1;		// Set error flag
    
				// Conservation of energy for internal
				// degrees of freedom
    dE             = 0.5*Mu*(*cr)*(*cr);

				// Remaining energy split between
				// internal degrees of freedom
    if (use_Eint>=0) {
      double kEm     = (totE - kE + dE)/(ne1 + ne2);
      p1->dattrib[use_Eint] = ne1 * kEm;
      p2->dattrib[use_Eint] = ne2 * kEm;
    }
    
				// Reset internal energy excess
    if (use_exes>=0) {
      if (ENSEXES) 		// Energy will be spread later
	p1->dattrib[use_exes] = p2->dattrib[use_exes] = 0.0;
      else {			// Energy excess incorporated now
	p1->dattrib[use_exes] =  p1->mass*(remE - delE)/Mt;
	p2->dattrib[use_exes] =  p2->mass*(remE - delE)/Mt;
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


void CollideIon::printCollInitialize() 
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
	    << "# ffN           number of free-free      " << std::endl
	    << "# ffE           cum energy in free-free  " << std::endl
	    << "# CEN           number of collions       " << std::endl
	    << "# CEE           cum energy in collisions " << std::endl
	    << "# CIN           number of ionizations    " << std::endl
	    << "# CIE           cum energy in ionizations" << std::endl
	    << "# RRN           number of rad recombs    " << std::endl
	    << "# RRE           energy in rad recombs    " << std::endl
	    << "# dEx           mean energy excess       " << std::endl
	    << "#"                                         << std::endl;

	out << "#" << std::right
	    << std::setw(11) << "Time "
	    << std::setw(12) << "ffN "
	    << std::setw(12) << "ffE "
	    << std::setw(12) << "CEN "
	    << std::setw(12) << "CEE "
	    << std::setw(12) << "CIN "
	    << std::setw(12) << "CIE "
	    << std::setw(12) << "RRN "
	    << std::setw(12) << "RRE "
	    << std::setw(12) << "dEx "
	    << std::endl;
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
	    << "# av            mean thermal energy      " << std::endl
	    << "# min           minimum thermal energy   " << std::endl
	    << "# max           maximum thermal energy   " << std::endl
	    << "# over10        number over 10.2 eV      " << std::endl
	    << "#"                                         << std::endl;

	out << "#" << std::right
	    << std::setw(11) << "Time "
	    << std::setw(12) << "av "
	    << std::setw(12) << "min "
	    << std::setw(12) << "max "
	    << std::setw(12) << "over10 "
	    << std::endl;
	}
      }
      in.close();
  }

}
   


void CollideIon::printCollSummary() 
{
  {
    std::ofstream out(coll_file_debug.c_str(), ios::out | ios::app);
    if (out) {
      out << std::setw(12) << tnow
	  << std::setw(12) << ff_d.first
	  << std::setw(12) << ff_d.second
	  << std::setw(12) << CE_d.first
	  << std::setw(12) << CE_d.second
	  << std::setw(12) << CI_d.first
	  << std::setw(12) << CI_d.second
	  << std::setw(12) << RR_d.first
	  << std::setw(12) << RR_d.second;
      if (dv.first>0.0)
	out << std::setw(12) << dv.second/dv.first << std::endl;
      else
	out << std::setw(12) << 0.0                << std::endl;
    }
  }

  {
    std::ofstream out(energy_file_debug.c_str(), ios::out | ios::app);
    if (out) {
      out << std::setw(12) << tnow
	  << std::setw(12) << eV_av/eV_N
	  << std::setw(12) << eV_min
	  << std::setw(12) << eV_max
	  << std::setw(12) << eV_10
	  << std::endl;
    }
  }
}

void CollideIon::resetColls() 
{
  ff_d    = std::make_pair(0., 0.);
  CE_d    = std::make_pair(0., 0.);
  CI_d    = std::make_pair(0., 0.);
  RR_d    = std::make_pair(0., 0.);
  dv      = std::make_pair(0., 0.);
  
  eV_N    = 0.0; 
  eV_av   = 0.0; 
  eV_max  = 0.0; 
  eV_min  = 99999.0; 

  eV_10   = 0.0;
}

void CollideIon::printCollGather() 
{
  if (myid==0) {
    MPI_Reduce(MPI_IN_PLACE, &ff_d.first,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 1
    MPI_Reduce(MPI_IN_PLACE, &ff_d.second, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 2
    MPI_Reduce(MPI_IN_PLACE, &CE_d.first,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 3
    MPI_Reduce(MPI_IN_PLACE, &CE_d.second, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 4
    MPI_Reduce(MPI_IN_PLACE, &CI_d.first,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 5
    MPI_Reduce(MPI_IN_PLACE, &CI_d.second, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 6
    MPI_Reduce(MPI_IN_PLACE, &RR_d.first,  1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 7
    MPI_Reduce(MPI_IN_PLACE, &RR_d.second, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 8
    MPI_Reduce(MPI_IN_PLACE, &dv.first,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 9
    MPI_Reduce(MPI_IN_PLACE, &dv.second,   1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 10
    MPI_Reduce(MPI_IN_PLACE, &eV_av,       1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 11
    MPI_Reduce(MPI_IN_PLACE, &eV_N,        1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 12
    MPI_Reduce(MPI_IN_PLACE, &eV_min,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 13
    MPI_Reduce(MPI_IN_PLACE, &eV_max,      1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 14
    MPI_Reduce(MPI_IN_PLACE, &eV_10,       1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 15
  } else {
    MPI_Reduce(&ff_d.first,  0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 1
    MPI_Reduce(&ff_d.second, 0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 2
    MPI_Reduce(&CE_d.first,  0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 3
    MPI_Reduce(&CE_d.second, 0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 4
    MPI_Reduce(&CI_d.first,  0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 5
    MPI_Reduce(&CI_d.second, 0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 6
    MPI_Reduce(&RR_d.first,  0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 7
    MPI_Reduce(&RR_d.second, 0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 8
    MPI_Reduce(&dv.first,    0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 9
    MPI_Reduce(&dv.second,   0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 10
    MPI_Reduce(&eV_av,       0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 11
    MPI_Reduce(&eV_N,        0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 12
    MPI_Reduce(&eV_min,      0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 13
    MPI_Reduce(&eV_max,      0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 14
    MPI_Reduce(&eV_10,       0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // 15
  }
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

