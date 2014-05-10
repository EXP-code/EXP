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

using namespace std;

// Proton mass (g)
const double mp      = 1.67262158e-24;

// Electron mass (g)
const double me      =  9.1094e-28;

// Speed of light (cm/s)
const double c       =  2.99792458e+10;

// electron volt in (cgs)
const double eV      =  1.6020e-12;

// Boltzmann constant (cgs)
const double boltz   = 1.3810e-16;

// Boltzmann constant (eV)
const double boltzEv = 8.6173324e-5;

// Planck's constant (cgs)
//const double planck  = 6.6262e-27;

// Electron charge (cgs)
const double esu     = 4.8032e-10;

//atomic mass unit in grams
const double amu = 1.660539e-24;

// Number of Ion elements
const int N_Z = 2;		

// Array of the different Z's for the elements
const int ZList[N_Z] = {1, 2}; 

// Harcoded abundance array
const double abund[N_Z] = {0.97, 0.03};

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
void CollideIon::initialize_cell(pHOT* tree, pCell* cell,
				 double rvmax, double tau,
				 sKey2Umap& nsel, int id)
{
  double KEtot, KEdspC;
  cell->KE(KEtot, KEdspC);	// KE in cell

  double massC = cell->Mass();	// Mass in cell

				// Used for diagnostics only
  totalSoFar += massC * KEdspC;
  massSoFar  += massC;
				// Representative avg cell velocity in cgs
  double vavg = 0.5*rvmax*UserTreeDSMC::Vunit;

				// Representative avg cell energy in ergs
  double Eerg = 0.5*vavg*vavg*amu/eV;

				// Upscaling factor for scattering cross section
  double sUp  = diamfac*diamfac;

  for (sKey2Umap::iterator it1 = nsel.begin(); it1 != nsel.end(); it1++)  {

    speciesKey i1 = it1->first;
    double Cross1 = geometric(i1.first);

    for (sKeyUmap::iterator 
	   it2 = it1->second.begin(); it2 != it1->second.end(); it2++)  
      {
	speciesKey i2 = it2->first;

	double Cross2 = geometric(i2.first);

	if (i1.second==1 && i2.second>1) 
	  Cross2 += elastic(i1.first, Eerg * atomic_weights[i1.first]) *
	    1e-14 / (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);

	csections[id][i1][i2] = (Cross1 + Cross2) * sUp;
      }
  }
}

sKey2Dmap& CollideIon::totalCrossSections(double crm, int id)
{    
  return csections[id];
}

sKey2Dmap& CollideIon::totalScatteringCrossSections(double crm, pCell *c, int id)
{
  typedef std::map<speciesKey, unsigned> Count;
  Count::iterator it1, it2;
  
  double vel = crm * UserTreeDSMC::Vunit;
  double Eerg = 0.5*vel*vel*amu/eV;

  // Upscaling factor for scattering cross section
  double sUp  = diamfac*diamfac;

  for (it1 = c->count.begin(); it1 != c->count.end(); it1++)  {

    speciesKey i1 = it1->first;
    double Cross1 = geometric(i1.first);
    
    for (it2 = c->count.begin(); it2 != c->count.end(); it2++)  
      {
	speciesKey i2 = it2->first;

	double Cross2 = geometric(i2.first);

	if (i1.second==1 && i2.second>1)
	  Cross2 += elastic(i1.first, Eerg * atomic_weights[i1.first]) *
	    1e-14 / (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
	
	csections[id][i1][i2] = (Cross1 + Cross2) * sUp;
      }
  }
    
  return csections[id];
}

int CollideIon::inelastic(pHOT *tree, Particle* p1, Particle* p2, 
			  double *cr, int id)
{
  int ret = 0;			// No error (flag)
  int interFlag = -1;

  // Number of atoms in each super particle
  //
  double N1 = (p1->mass*UserTreeDSMC::Munit)/(atomic_weights[p1->Z]*amu);
  double N2 = (p2->mass*UserTreeDSMC::Munit)/(atomic_weights[p2->Z]*amu);
  
  // Number of associated electrons for each particle scaled by the
  // number in the atomic target
  //
  double ne1 = N1/N2 * (p1->C - 1);
  double ne2 = N2/N1 * (p2->C - 1);
  //           ^       ^
  //           |       |
  //           |       +---Number of FREE electrons per atom
  //           |
  // Statistical weight owing to number of particles per superparticle
  //

  // The total mass in system units
  //
  double Mt = p1->mass + p2->mass;
  if (Mt<=0.0) return ret;
  
  // Reduced mass in ballistic collision (system)
  //
  double Mu = p1->mass * p2->mass /Mt;

  // Center of mass energy in the ballistic collision (system)
  double kE  = 0.5*Mu*(*cr)*(*cr);

  // Energy available in the center of mass of the atomic collision
  //
  double m1  = atomic_weights[p1->Z]*amu;
  double m2  = atomic_weights[p2->Z]*amu;
  double mu  = m1 * m2 / (m1 + m2);
  double kEI = 0.5*mu*(*cr)*(*cr) * UserTreeDSMC::Vunit*UserTreeDSMC::Vunit;

  // Convert ergs to eV
  //
  double kEe = kEI * 6.241509e11; 
  
  //
  // For tracking energy consistency
  //
  double dE   = kE*TolV*TolV;
  double remE = (kE - dE);
  double delE = 0.0, delEeV = 0.0;

  
  /** std::cout << "total Ke = " << kE << " with cr = " << (*cr) 
      << " dE = " << dE << " remE = " << remE << endl; */
  
  // Get temperatures from cells
  //
  double E_therm_1 = boltzEv*p1->dattrib[use_temp];
  double E_therm_2 = boltzEv*p2->dattrib[use_temp];
  
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
  std::vector<double> dCross;

  // Index the interactions
  std::vector<int> inter;

  double sum12 = 0.0;		// Accumulate total cross sections as we go
  double sum21 = 0.0;
  
  // Energy summary diagnostics
  eV_av += kEe;
  if (std::isnan(eV_av)) {
    std::cout << "eV_N=" << eV_N << std::endl;
  }
  eV_N++;
  eV_min = std::min(eV_min, kEe);
  eV_max = std::max(eV_max, kEe);

  if (kEe > 10.2) { eV_10++;}
  
  //--------------------------------------------------
  // Particle 1 interacts with Particle 2
  //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
  if (p1->C > 1 and ne2 > 0) {
    double ff1 = IonList[p1->Z][p1->C].freeFreeCross(ch, kEe);

    dCross.push_back(ne2*ff1);
    sum12 += ff1*ne2;
    inter.push_back(1);
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  std::vector< std::pair<double, double > > CE1;
  if (ne2 > 0 and p1->C <= p1->Z) {

    assert(E_therm_1 != 0);
    
    CE1 = IonList[p1->Z][p1->C].collExciteCross(ch, kEe, E_therm_1);

    dCross.push_back(ne2*CE1.back().first);
    sum12 += CE1.back().first*ne2;
    inter.push_back(2);
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  double DI1;
  if (p1->C < (p1->Z + 1) and ne2 > 0) {

    DI1 = IonList[p1->Z][p1->C].directIonCross(ch, kEe);

    dCross.push_back(ne2*DI1);
    sum12 += DI1*ne2;
    inter.push_back(3);
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  std::vector<double> RE1;
  if (p1->C > 1 and ne2 > 0) {

    RE1 = IonList[p1->Z][p1->C].radRecombCross(ch, kEe);

    dCross.push_back(ne2*RE1.back());
    sum12 += RE1.back()*ne2;
    inter.push_back(4);
  }

  
  //--------------------------------------------------
  // Particle 2 interacts with Particle 1
  //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
  if (p2->C > 1 and ne1 > 0) {
    double ff2 = IonList[p2->Z][p2->C].freeFreeCross(ch, kEe);

    dCross.push_back(ne1*ff2);
    sum21 += ff2*ne1;
    inter.push_back(6);
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  std::vector< std::pair<double, double > > CE2;
  if(ne1 > 0 and p2->C <= p2->Z) {
    assert(E_therm_2 != 0);
    CE2 = IonList[p2->Z][p2->C].collExciteCross(ch, kEe, E_therm_2);

    dCross.push_back(ne1*CE2.back().first);
    sum21 += CE2.back().first*ne1;
    inter.push_back(7);
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  double DI2;
  if (p2->C < (p2->Z + 1) and ne1 > 0) {
    DI2 = IonList[p2->Z][p2->C].directIonCross(ch, kEe);

    dCross.push_back(ne1*DI2);
    sum21 += DI2*ne1;
    inter.push_back(8);
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  std::vector<double> RE2;
  if (p2->C > 1 and ne1 > 0) {
    RE2 = IonList[p2->Z][p2->C].radRecombCross(ch, kEe);
    dCross.push_back(ne1*RE2.back());

    sum21 += RE2.back()*ne1;
    inter.push_back(9);
  } 
  
  // Now that the interactions have been calculated, create the
  // normalized cross section list to pick the interaction

  std::vector<double> TotalCross;
  double tCross;
  tCross = 0;
  int si = 0;
  for(size_t i = 0; i < dCross.size(); i++) {
    tCross += dCross[i];
    TotalCross.push_back(tCross);
    si++;
  }
  assert(TotalCross.size() == dCross.size());

  if (tCross != 0) {
    std::vector<double> CDF;
    std::vector <double>::iterator 
      max = max_element(TotalCross.begin(), TotalCross.end());

    for (size_t i = 0; i < TotalCross.size(); i++) {
      CDF.push_back(TotalCross[i]/tCross);
    }
    
    // Use a random variate to select the interaction from the
    // discrete cumulatative probability distribution (CDF)
    double ran = (*unit)();
    int index  = -1;
    for (size_t i = 0; i < CDF.size(); i++) {
      if (ran < CDF[i]) {
	index = static_cast<int>(i);
	break;
      }
    }

    // Sanity check: did not assign index??
    if (index<0) {
      std::cout << "CDF location falure, myid=" << myid
		<< ", ran=" << ran << std::endl;
      index = 0;
    }

    // Finally, set the interaction type based on the selected index
    interFlag = inter[index];

    // Sanity check: did not set interaction type??
    if (interFlag<0) {
      std::cout << "interFlag NOT set, myid=" << myid 
		<< ", index=" << index << std::endl;
      index = 0;
      interFlag = inter[0];
    }

    int partflag = 0;		// Will be 1 or 2, dependending on
				// which ion or neutral is selected
				// for interaction

    //-------------------------
    // Particle 1 interactions
    //-------------------------

    if (interFlag == 1) {
      delE          = IS.selectFFInteract(IonList[p1->Z][p1->C], kEe);
      partflag      = 1;
      ff_d.first++; 
      ff_d.second  += delE;
    }

    if (interFlag == 2) {
      delE = IS.selectCEInteract(IonList[p1->Z][p1->C], CE1);
      partflag      = 1;
      CE_d.first++; 
      CE_d.second  += delE;
    }

    if (interFlag == 3) {
      delE          = IS.DIInterLoss(ch, IonList[p1->Z][p1->C]);
      p1->C++;
      assert(p1->C <= (p1->Z + 1));
      partflag      = 1;
      CI_d.first++; 
      CI_d.second  += delE;
    }

    // KE carried by electron is subtracted from the thermal reservoir
    // but not the radiation of "binding".  This radiation decreases
    // the total energy of the gas but not the thermal component.
    if (interFlag == 4) {
      delE          = kEe;
      p1->C--;
      assert(p1->C > 0);
      partflag      = 1;
      RR_d.first++; 
      RR_d.second  += delE;
    }
    
    //-------------------------
    // Particle 2 interactions
    //-------------------------

    if (interFlag == 6) {
      delE          = IS.selectFFInteract(IonList[p2->Z][p2->C], kEe);
      partflag      = 2;
      ff_d.first++;
      ff_d.second += delE;
    }

    if (interFlag == 7) {
      delE         = IS.selectCEInteract(IonList[p2->Z][p2->C], CE2);
      partflag     = 2;
      CE_d.first++; 
      CE_d.second += delE;
    }

    if (interFlag == 8) {
      delE = IS.DIInterLoss(ch, IonList[p2->Z][p2->C]);
      p2->C++;
      assert(p2->C <= (p2->Z + 1));
      CI_d.first++; 
      CI_d.second += delE;
      partflag     = 2;
    }

    if (interFlag == 9) {
      delE         = kEe;	// See comment above for interFlag==4
      p2->C--;
      assert(p2->C > 0);
      partflag     = 2;
      RR_d.first++; 
      RR_d.second += delE;
    }

    // Convert from eV to system units
    // Get the number of particles in cell and mult. by delE

    delEeV = delE;

    if (partflag == 1) delE *= N1;
    if (partflag == 2) delE *= N2;
    
    // Convert back to cgs
    //
    delE = delE*1.602177e-12;
  }
  
  assert(delE >= 0.);
  
  // concert to system units
  delE = delE/UserTreeDSMC::Eunit;
  
  if (NO_COOL) {
    delE = 0.;
  }
  
  if (frost_warning && delE*0.999 > kE) 
    std::cout << "delE > KE!! (" << delE << " > " << kE
	      << "), Interaction type = " << interFlag 
	      << " kEe = "  << kEe
	      << " delE = " << delEeV 
	      << std::endl;
  
  // Add the cross sections into the csections[id]
  // HERE
  speciesKey j1(p1->Z, p1->C);
  speciesKey j2(p2->Z, p2->C);
  
  // Add to the total elastic scattering cross section to the inelastic
  // contribution from the "interaction" section

  // Convert from nanometer^2 to cm^2 first, then system units

  double sUp = diamfac*diamfac;

  csections[id][j2][j1] = (sum12 + elastic(p1->Z, kEe)*sUp) * 1.0e-14 /
    (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);

  csections[id][j1][j2] = (sum21 + elastic(p2->Z, kEe)*sUp) * 1.0e-14 / 
    (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
  
  if (remE<=0.0 || delE<=0.0) return ret;
  
  // Cooling rate diagnostic

  if (TSDIAG) {
    if (delE>0.0) {
      int indx = (int)floor(log(remE/delE)/(log(2.0)*TSPOW) + 5);
      if (indx<0 ) indx = 0;
      if (indx>10) indx = 10;
      
      EoverT[id][indx] += Mt;
    }
  }
  
  if (use_exes>=0) {
    // (-/+) value means under/overcooled: 
    // positive/negative increment to delE
    // NB: delE may be < 0 if too much energy 
    // was radiated previously . . .
    //
    delE -= p1->dattrib[use_exes] + p2->dattrib[use_exes];
  }
  
  if (remE > delE) {
    double vi      = (*cr);
    lostSoFar[id] += delE;
    decelT[id]    += delE;
    (*cr)          = sqrt( 2.0*(kE - delE)/Mu );
    dv.first++; 
    dv.second     += 0.5*(vi - (*cr))*(vi - (*cr));
    ret            = 0;		// No error
    
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

