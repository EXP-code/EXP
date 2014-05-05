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

bool vdebug = false;

CollideIon::CollideIon(ExternalForce *force, double diameter, int Nth) : 
  Collide(force, diameter, Nth)
{
  NUM = 0;
  std::cout << "Creating CollideIon instance" << endl;
  csections = std::vector<sKey2Dmap> (nthrds);
  
  std::cout << "\tCreating ion list" << endl;
  for(int i = 0; i < N_Z; i++) {
    for (int j = 1; j <= ZList[i] + 1; j++) {
      IonList[ZList[i]][j] = Ion(ZList[i], j, ch);
      IonList[ZList[i]][j].freeFreeDifferential(ch);
    }
  }

  for (int i = 0; i < N_Z; i++) {
    Ni.push_back(1); // dummy value for now
  }
  firstFlag = 0;
  
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
  atomic_weights[1] = 1.0079;
  atomic_weights[2] = 4.0026;
  atomic_weights[3] = 6.941;
  atomic_weights[4] = 9.0122;
  atomic_weights[5] = 10.811;
  atomic_weights[6] = 12.011;
  atomic_weights[7] = 14.007;
  atomic_weights[8] = 15.999;
  atomic_weights[9] = 18.998;
  atomic_weights[10] = 20.180;
  atomic_weights[11] = 22.990;
  atomic_weights[12] = 24.305;
  
  ff_d   = std::make_pair(0.,0.);
  CE_d   = std::make_pair(0.,0.);
  CI_d   = std::make_pair(0.,0.);
  RR_d   = std::make_pair(0.,0.);
  dv     = std::make_pair(0., 0.);
  eV_av  = 0;
  eV_N   = 0;
  eV_min = 999999;
  eV_max = 0;
  eV_10  = 0;
  
  std::cout << "No cooling is set to: " << NO_COOL << endl;
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
  // std::cout << "Initializing cell in CollideIon " << endl;

  if (firstFlag == 0) {
    double cross = M_PI*diam*diam;
    for (sKey2Umap::iterator it1 = nsel.begin(); it1 != nsel.end(); it1++) 
      {
	speciesKey i1 = it1->first;
	for (sKeyUmap::iterator 
	       it2 = it1->second.begin(); it2 != it1->second.end(); it2++) 
	  {
	    speciesKey i2 = it2->first;
	    csections[id][i1][i2] = 
	      0.25*(i1.first*i1.first + i2.first*i2.first)*cross;
	  }
      }
    
    // firstFlag = 1;
  }
}

sKey2Dmap& CollideIon::totalCrossSections(double crm, int id)
{
  return csections[id];
}

int CollideIon::inelastic(pHOT *tree, Particle* p1, Particle* p2, 
			  double *cr, int id)
{
  int ret = 0;			// No error (flag)
  int interFlag = -1;

  // Number of associated electrons for each particle
  double ne1 = (p1->C - 1);
  double ne2 = (p2->C - 1);
  
  double N1 = (p1->mass*UserTreeDSMC::Munit)/(atomic_weights[p1->Z]*amu);
  double N2 = (p2->mass*UserTreeDSMC::Munit)/(atomic_weights[p2->Z]*amu);	
  
  // The relative velocity in physical units
  double vr = (*cr)*UserTreeDSMC::Vunit;
  
  // Physical masses for each particle
  double m1 = atomic_weights[p1->Z]*amu;
  double m2 = atomic_weights[p2->Z]*amu;
  
  // The reduced mass
  double muI = (m1*m2)/(m1 + m2);
  
  // The total mass in system units
  double Mt = p1->mass + p2->mass;
  if (Mt<=0.0) return ret;
  
  // The reduced mass in system units
  double Mu = p1->mass*p2->mass/Mt;
  if (Mu<=0.0) return ret;
  
  // Energy floor
  double kE  = 0.5*Mu*(*cr)*(*cr); // System units
  double kEI = 0.5*muI*vr*vr;	   // Physical units

  /**
     The physical approximation
     --------------------------
     The COM energy is shared, on average, among all ions and electrons
  */

  // Factor in the energy spread to all the particles
  kEI *= 2.0/(2.0+ne1+ne2);

  // The average electron velocity in this shared-energy scenario
  double pe = sqrt(2.*me*kEI);

  // std::cout << "Electron velocity = " << pe/me << endl;

  double kEe;
  kEe = kEI;
  kEe = kEe*6.241509e11;	// Convert ergs to eV
  
  //
  // For tracking energy consistency
  //
  double dE   = kE*TolV*TolV;
  double remE = (kE - dE);
  double delE = 0.0;
  
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
  
  // What is eV_av?
  eV_av += kEe;
  eV_N++;
  eV_max = max(eV_max, kEe);
  eV_min = min(eV_min, kEe);
  double delEeV;
  if (kEe > 10.2) { eV_10++;}
  
  int outflag = 1;
  
  // cout << "Calculating particle interactions in inelastic " << endl;
  // cout << "Cross-sections: " << endl;
  
  //--------------------------------------------------
  // Particle 1 interactions with 2
  //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
  if (p1->C > 1 and ne2 > 0) {
    double ff1 = IonList[p1->Z][p1->C].freeFreeCross(ch, kEe);
    dCross.push_back(ne2*ff1);
    // if (ff1 != 0) cout << "FF1: " << ff1 << " kE = " << kEe << endl;
    sum12 += ff1*ne2;
    inter.push_back(1);
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  std::vector< std::pair<double, double > > CE1;
  if (ne2 > 0 and p1->C <= p1->Z) {

    assert(E_therm_1 != 0);
    
    CE1 = IonList[p1->Z][p1->C].collExciteCross(ch, pe, E_therm_1, me);

    /**
       if (ne2*CE1.back().first != 0) {
          std::cout << "\tCE1: " << CE1.back().first*ne2 << endl;
       }
    */

    /** if (CE1.back().first != 0) 
	   cout << "CE = " << ne2*CE1.back().first << " E = " << E_therm_1 
	        << " Z = " << p1->Z << " C = " << p1->C <<"\t";
    */
    
    dCross.push_back(ne2*CE1.back().first);
    sum12 += CE1.back().first*ne2;
    inter.push_back(2);
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  double DI1;
  // std::cout << "THREE\t";
  if (p1->C < (p1->Z + 1) and ne2 > 0) {
    DI1 = IonList[p1->Z][p1->C].directIonCross(ch, kEe);
    /**
     if (DI1 != 0) { 
       cout << "\tDI1 = " << DI1 << endl;
    }
    */
    dCross.push_back(ne2*DI1);
    sum12 += DI1*ne2;
    inter.push_back(3);
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  std::vector<double> RE1;
  // std::cout << "FOURTH\t";
  if (p1->C > 1 and ne2 > 0) {
    RE1 = IonList[p1->Z][p1->C].radRecombCross(ch, kEe);
    dCross.push_back(ne2*RE1.back());
    outflag = 0;
    // cout << "RR1 = " << ne2*RE1.back() <<  " kE = " << kEe << endl;
    sum12 += RE1.back()*ne2;
    inter.push_back(4);
  }

  
  //--------------------------------------------------
  // Particle 2 interactions with 1
  //--------------------------------------------------

				//-------------------------------
				// *** Free-free
				//-------------------------------
  // std::cout << "SIXTH\t";
  if (p2->C > 1 and ne1 > 0) {
    double ff2 = IonList[p2->Z][p2->C].freeFreeCross(ch, kEe);
    /**
       if (ff2 != 0) 
         cout << "FF2: " << ff2 << " kE = " << kEe << endl;
    */
    dCross.push_back(ne1*ff2);
    sum21 += ff2*ne1;
    inter.push_back(6);
  }
				//-------------------------------
				// *** Collisional excitation
				//-------------------------------
  // std::cout << "SEVENTH\t (Z,C) = " << p2->Z << "," << p2->C << ")" << "\t";
  std::vector< std::pair<double, double > > CE2;
  if(ne1 > 0 and p2->C <= p2->Z) {
    assert(E_therm_2 != 0);
    CE2 = IonList[p2->Z][p2->C].collExciteCross(ch, pe, E_therm_2, me);
    /**
       if (ne1*CE2.back().first != 0) {
          cout << "\tCE2: " << CE2.back().first*ne1 << endl;
       }
       if (CE2.back().first != 0) 
          cout << "CE = " << ne1*CE2.back().first << "\t";
    */
    dCross.push_back(ne1*CE2.back().first);
    sum21 += CE2.back().first*ne1;
    inter.push_back(7);
  }
				//-------------------------------
				// *** Ionization cross section
				//-------------------------------
  double DI2;
  // std::cout << "EIGHTH\t";
  if (p2->C < (p2->Z + 1) and ne1 > 0) {
    DI2 = IonList[p2->Z][p2->C].directIonCross(ch, kEe);
    /**
       if (DI2 != 0) {
          cout << "\tDI2 = " << DI2 << endl;
       }
    */
    dCross.push_back(ne1*DI2);
    sum21 += DI2*ne1;
    inter.push_back(8);
  }
				//-------------------------------
				// *** Radiative recombination
				//-------------------------------
  std::vector<double> RE2;
  // std::cout << "NINETH\n";
  if (p2->C > 1 and ne1 > 0) {
    RE2 = IonList[p2->Z][p2->C].radRecombCross(ch, kEe);
    dCross.push_back(ne1*RE2.back());
    // std::cout << "RR2 = " << ne1*RE2.back() << " kE = " << kEe << endl;
    sum21 += RE2.back()*ne1;
    inter.push_back(9);
    outflag = 0;
  } 
  
  // std::cout << "DONE\t";

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
    std::vector<double> normed;
    std::vector <double>::iterator 
      max = max_element(TotalCross.begin(), TotalCross.end());

    for (size_t i = 0; i < TotalCross.size(); i++) {
      normed.push_back(TotalCross[i]/tCross);
    }
    
    // Now that the normed probablility cross section vector is
    // created, pull a random number
    double ran = (*unit)();
    int index  = -1;

    // Locate the interaction in the cumulative cross-section list
    for (size_t i = 0; i < normed.size(); i++) {
      if (ran < normed[i]) {
	index = static_cast<int>(i);
	break;
      }
    }

    // Sanity check
    if (index<0) {
      std::cout << "CDF location falure, ran=" << ran << std::endl;
      index = 0;
    }

    if (vdebug) {
      if (outflag) std::cout << index << "\t";
      interFlag = inter[index];
    }

    // Sanity check
    if (interFlag<0) {
      std::cout << "interFlag NOT set, index=" << index << std::endl;
      interFlag = inter[0];
    }

    int partflag = 0;

    if (vdebug) {
      if (outflag) std::cout << interFlag << std::endl;
      std::cout << "Interactions\n";
    }

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
      if (vdebug) cout << "\nExcitation 1" << endl;
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

    if (interFlag == 4) {
      delE          = 0;
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
      if (vdebug) cout << "\nExcitation 2" << endl;
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
      delE         = 0.0;
      p2->C--;
      assert(p2->C > 0);
      partflag     = 2;
      RR_d.first++; 
      RR_d.second += delE;
    }

    // Convert from eV to system units
    // Get the number of particles in cell and mult. by delE

    double NP;
    if (vdebug) {
      if (partflag == 1) 
	NP = N1;
      if (partflag == 2)
	NP = N2;
    }
    
    NP = Mu*UserTreeDSMC::Munit/muI;
    delEeV = delE;
    delE   = delE*NP;
    
    // convert to cgs
    delE = delE*1.602177e-12;
    
  }
  
  assert(delE >= 0.);
  
  // concert to system units
  delE = delE/UserTreeDSMC::Eunit;
  
  if(NO_COOL) {
    delE = 0.;
  }
  
  if (delE > kE) 
    std::cout << "delE > KE!! Interaction = " << interFlag 
	      << " kEe = " << kEe << " kE = " << kE 
	      << " delE = " << delEeV << std::endl;
  /**
  if (delE != 0.) 
    std::cout << "delE = " << delE << "kE = " << kE 
	      << " remE = " << remE << std::endl;
  */
  
  
  // Add the cross sections into the csections[id]
  // HERE
  speciesKey j1(p1->Z, p1->C);
  speciesKey j2(p2->Z, p2->C);
  
  // add to the total cross section the "ballistic" minimum bohr radius
  // for the cross section. THe first part is then the "interaction"
  // section, and the second is the elastic collision

  // convert to cm^2 first, then system units

  csections[id][j2][j1] = 
    (sum12*1e-14)/(UserTreeDSMC::Lunit*UserTreeDSMC::Lunit) + 
    0.25*(p1->Z*p1->Z+p2->Z*p2->Z)*M_PI*diam*diam; 

  csections[id][j1][j2] = 
    (sum21*1e-14)/(UserTreeDSMC::Lunit*UserTreeDSMC::Lunit) +
    0.25*(p1->Z*p1->Z+p2->Z*p2->Z)*M_PI*diam*diam;
  
  // std::cout << "dE = " << delE;
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
  
  if (remE >= delE) {
    double vi      = (*cr);
    lostSoFar[id] += delE;	
    decelT[id]    += delE;
    (*cr)          = sqrt( 2.0*(kE - delE)/Mu );
    dv.first++; 
    dv.second      = vi - (*cr);
    ret            = 0;		// No error
    
    /*
    if (interFlag == 1 or interFlag == 6)
      std::cout << "Before: " << vi 
		<< " After: " << (*cr) 
		<< " dE = " << delE << endl;
    */

				// Zero out internal energy excess
    if (use_exes>=0)		// since excess is now used up
      p1->dattrib[use_exes] = p2->dattrib[use_exes] = 0.0;

  } else {			// Inconsistent: too much energy lost!

    std::cout << "************remE < delE!!!!" 
	      << " RemE = " << remE << " delE = " << delE 
	      << " interaction: " << interFlag    << std::endl;
    
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

void CollideIon::printCollSummary() 
{
  if (myid == 0) {
    std::cout << "Collision debugging:"                  << std::endl;
    std::cout << "id:ffN:ffE:CEN:CEE:CIN:CIE:RRN:RRE:dv" << std::endl;
    std::cout << myid  << ff_d.first << " : " << ff_d.second 
	      << " : " << CE_d.first << " : " << CE_d.second 
	      << " : " << CI_d.first << " : " << CI_d.second 
	      << " : " << RR_d.first << " : " << RR_d.second 
	      << " : " << dv.second/dv.first  << std::endl;

  } else { 
    std::cout << myid  << ff_d.first << " : " << ff_d.second 
	      << " : " << CE_d.first << " : " << CE_d.second 
	      << " : " << CI_d.first << " : " << CI_d.second 
	      << " : " << RR_d.first << " : " << RR_d.second 
	      << " : " << dv.second/dv.first  << std::endl; 
  }
  
  if (myid == 0) {
    std::cout << "Energy debugging:"    << std::endl;
    std::cout << "id:av:min:max:over10" << std::endl;
    std::cout << myid 
	 << " : " << eV_av/eV_N << " : " << eV_min 
	 << " : " << eV_max     << " : " << eV_10 << std::endl;
  }
  else { 
    std::cout << myid 
	      << " : " << eV_av/eV_N  
	      << " : " << eV_min 
	      << " : " << eV_max << endl; 
  }
}

void CollideIon::resetColls() 
{
  ff_d    = std::make_pair(0.,0.);
  CE_d    = std::make_pair(0.,0.);
  CI_d    = std::make_pair(0.,0.);
  RR_d    = std::make_pair(0.,0.);
  dv      = std::make_pair(0., 0.);
  
  eV_N    = 0.0; 
  eV_av   = 0.0; 
  eV_max  = 0.0; 
  eV_min  = 99999.0; 

  eV_10++;
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
    
    crossIJ = totalCrossSections(0, id);
    
    for (it1=c->count.begin(); it1!=c->count.end(); it1++) {
      speciesKey i1 = it1->first;
      densM[i1] = c->Mass(i1)/volc;
    }
    
    double meanDens=0.0, meanLambda=0.0;
    
    for (it1=c->count.begin(); it1!=c->count.end(); it1++) {
      speciesKey i1 = it1->first;
      /*
	std::cout << "species: " 
	<< "( " << i1.first << " , " << i1.second << " ) " << std::endl;
      */
      crossM [i1] = 0.0;
      for (it2=c->count.begin(); it2!=c->count.end(); it2++) {
	speciesKey i2 = it2->first;
	double      N = UserTreeDSMC::Munit/(amu*atomic_weights[i2.first]);
	
	/*
	std::cout << "N = " << N << " dens = " << densM[i2] << std::endl;
	*/
	crossIJ[i1][i2] = std::max<double>
	  (
	   crossIJ[i1][i2], 
	   (M_PI*0.25*diam*diam*(i1.first*i1.first + i2.first*i2.first))
	   );

	if (i2>=i1) {
	  crossM[i1] += N*densM[i2]*std::min<double>
	    (
	     crossIJ[i1][i2], 
	     (M_PI*0.25*diam*diam*(i1.first*i1.first + i2.first*i2.first))
	     );

	  /*
	    std::cout << i1.first << " " << i1.second << " " 
	              << i2.first << " " << i2.second << " " 
		      << crossIJ[i1][i2] << std::endl;
	  */
	} else
	  crossM[i1] += N*densM[i2]*std::min<double>
	    (
	     crossIJ[i2][i1], 
	     (M_PI*0.25*diam*diam*(i1.first*i1.first + i2.first*i2.first))
	     );
	/*
	  std::cout << i2.first << " " << i2.second << " " 
	            << i1.first << " " << i1.second << " " 
		    << crossIJ[i2][i1] << std::endl;
	*/
      }
      
      /*
	std::cout << "crossM = " << crossM[i1] 
	          << " densM = " << densM[i1] 
		  << "Fn: " << (*Fn)[i1] << std::endl;
      */
      lambdaM[i1] = 1.0/crossM[i1];
      meanDens   += densM[i1] ;
      meanLambda += densM[i1] * lambdaM[i1];
    }
    
    // This is the number density-weighted
    meanLambda /= meanDens;
    
    /*
      std::cout << meanLambda << "\t" 
                << meanLambda*UserTreeDSMC::Lunit << std::endl;
    */
    
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

