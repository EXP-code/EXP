#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cfloat>
#include <cmath>
#include <ctime>
#include <tuple>
#include <map>

#include <boost/filesystem.hpp>

#include "global.H"
#include "UserTreeDSMC.H"
#include "CollideIon.H"
#include "localmpi.h"
#include "Species.H"
#include "Configuration.H"

using namespace std;

double   CollideIon::Nmin    = 1.0e-08;	   
double   CollideIon::Nmax    = 1.0e+25;	   
double   CollideIon::Tmin    = 1.0e+03;	   
double   CollideIon::Tmax    = 1.0e+08;	   
unsigned CollideIon::Nnum    = 400;	   
unsigned CollideIon::Tnum    = 200;	   
string   CollideIon::cache   = ".HeatCool";
bool     CollideIon::equiptn = false;
bool     CollideIon::scatter = false;
bool     CollideIon::ExactE  = false;
bool     CollideIon::AlgOrth = false;
bool     CollideIon::DebugE  = false;
unsigned CollideIon::esNum   = 100;
unsigned CollideIon::NoDelC  = 0;
double   CollideIon::logL    = 10;
string   CollideIon::config0 = "CollideIon.config";

CollideIon::ElectronScatter
CollideIon::esType           = CollideIon::always;

CollideIon::esMapType CollideIon::esMap = { {"none",      none},
					    {"always",    always},
					    {"classical", classical},
					    {"limited",   limited},
					    {"fixed",     fixed} };

// Add trace energy excess to electron distribution
//
static bool TRACE_ELEC        = false;

// Enable ion-electron secondary scattering
//
static bool SECONDARY_SCATTER = false;

// Fraction of excess energy loss to give to the electrons
//
static double TRACE_FRAC      = 1.0;

// Apply "splitting" energy to the COM scattering regardless of trace type
//
static bool TRACE_REAPPLY     = false;

// Print collisions by species for debugging
//
static bool COLL_SPECIES      = false;

// Same species tests (for debugging only)
//
static bool SAME_ELEC_SCAT    = false;
static bool SAME_IONS_SCAT    = false;
static bool SAME_INTERACT     = false;
static bool SAME_TRACE_SUPP   = false;

// Suppress distribution of energy to electrons when using NOCOOL
//
static bool NOCOOL_ELEC       = false;

// Suppress distribution of ionization energy between electrons
//
static bool NOSHARE_ELEC      = false;

// Clone temperature of ionizing electron
//
static bool CLONE_ELEC        = false;

// Warn if energy lost is smaller than COM energy available.  For
// debugging.  Set to false for production.
//
static bool frost_warning     = false;

// Very verbose selection debugging. Set to false for production.
//
static bool DEBUG_SL          = false;

// Verbose cross-section debugging. Set to false for production.
//
static bool DEBUG_CR          = false;

// Verbose cross-section debugging for unequal species only. Set to
// false for production.
//
static bool DEBUG_NQ          = false;

// Artifically suppress electron equipartition speed
//
static bool NO_DOF            = true;

// Artifically suppress electron equilibrium velocity
//
static bool NO_VEL            = false;

// Artifically suppress energy loss due to ionization
//
static bool NO_ION_E          = false;

// Artifically suppress energy loss due to free-free
//
static bool NO_FF_E          = false;

// KE debugging: checks energy bookkeeping for weighted algorithm. Set
// to false for production
//
static bool KE_DEBUG          = true;

// KE debugging threshold for triggering diagnostic output
//
static double DEBUG_THRESH    = 1.0e-9;

// Tally ionization potential with energy loss during recombination
//
static bool RECOMB_IP         = false;

// Cross-section debugging; set to false for production
//
static bool CROSS_DBG         = false;

// Excess trace map debugging; set to false for production
//
static bool EXCESS_DBG        = false;

// Enable NTC full distribution for electrons
static bool NTC_DIST          = true;

// Minimum energy for Rutherford scattering of ions used to estimate
// the elastic scattering cross section
//
static double FloorEv         = 0.05;

// Minimum relative fraction for allowing a collisional excitation
//
static double minCollFrac     = -1.0;

// Test use_cons summation for debugging
//
static bool use_cons_test     = true;

CollideIon::CollideIon(ExternalForce *force, Component *comp, 
		       double hD, double sD, 
		       const std::string& smap, int Nth) : 
  Collide(force, comp, hD, sD, Nth)
{
  // Process the feature config file
  //
  processConfig();

  // Debugging
  //
  itp=0;

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

  // Banners logging the test algorithms
  //
  if (myid==0 && NOCOOL)
    std::cout << std::endl
	      << "************************************" << std::endl
	      << "*** No cooling is ON for testing ***" << std::endl
	      << "************************************" << std::endl;

  if (myid==0 && scatter)
    std::cout << std::endl
	      << "************************************" << std::endl
	      << "*** No recombination/ionization  ***" << std::endl
	      << "************************************" << std::endl;

  if (myid==0 && equiptn)
    std::cout << std::endl
	      << "************************************" << std::endl
	      << "*** Using electron EQUIPARTITION ***" << std::endl
	      << "************************************" << std::endl;


  if (myid==0) {
    std::cout << std::endl
	      << "************************************" << std::endl
	      << "*** Algorithm selection flags ******" << std::endl
	      << "************************************" << std::endl
	      << " " << std::setw(20) << std::left  << "ENERGY_ES"
	      << (ExactE ? "on" : "off")                << std::endl
	      <<  " " << std::setw(20) << std::left << "ENERGY_DBG"
	      << (DebugE ? "on" : "off")                << std::endl
	      <<  " " << std::setw(20) << std::left << "ENERGY_ORTHO"
	      << (AlgOrth ? "on" : "off")               << std::endl
	      <<  " " << std::setw(20) << std::left << "SECONDARY_SCATTER"
	      << (SECONDARY_SCATTER ? "on" : "off")     << std::endl
	      <<  " " << std::setw(20) << std::left << "COLL_SPECIES"
	      << (COLL_SPECIES ? "on" : "off")         << std::endl
	      <<  " " << std::setw(20) << std::left << "TRACE_ELEC"
	      << (TRACE_ELEC ? "on" : "off")            << std::endl
	      <<  " " << std::setw(20) << std::left << "TRACE_FRAC"
	      << TRACE_FRAC                             << std::endl
	      <<  " " << std::setw(20) << std::left << "TRACE_REAPPLY"
	      << (TRACE_REAPPLY ? "on" : "off")         << std::endl
	      <<  " " << std::setw(20) << std::left << "SAME_ELEC_SCAT"
	      << (SAME_ELEC_SCAT ? "on" : "off")        << std::endl
	      <<  " " << std::setw(20) << std::left << "SAME_IONS_SCAT"
	      << (SAME_IONS_SCAT ? "on" : "off")        << std::endl
	      <<  " " << std::setw(20) << std::left << "SAME_INTERACT"
	      << (SAME_INTERACT ? "on" : "off")         << std::endl
	      <<  " " << std::setw(20) << std::left << "SAME_TRACE_SUPP"
	      << (SAME_TRACE_SUPP ? "on" : "off")       << std::endl
	      <<  " " << std::setw(20) << std::left << "NoDelC"
	      << NoDelC                                 << std::endl
	      <<  " " << std::setw(20) << std::left << "NOCOOL_ELEC"
	      << (NOCOOL_ELEC ? "on" : "off")           << std::endl
	      <<  " " << std::setw(20) << std::left << "NOSHARE_ELEC"
	      << (NOSHARE_ELEC ? "on" : "off")          << std::endl
	      <<  " " << std::setw(20) << std::left << "CLONE_ELEC"
	      << (CLONE_ELEC ? "on" : "off")            << std::endl
	      <<  " " << std::setw(20) << std::left << "RECOMB_KE"
	      << (RECOMB_IP ? "on" : "off")             << std::endl
	      <<  " " << std::setw(20) << std::left << "KE_DEBUG"
	      << (KE_DEBUG ? "on" : "off" )             << std::endl
	      <<  " " << std::setw(20) << std::left << "NTC_DIST"
	      << (NTC_DIST ? "on" : "off" )             << std::endl
	      << "************************************" << std::endl;
  }

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
  kEee     .resize(nthrds);
  Ein1     .resize(nthrds);
  Ein2     .resize(nthrds);
  Evel     .resize(nthrds);
  spTau    .resize(nthrds);
  spCrm    .resize(nthrds);
  spNsel   .resize(nthrds);
  spProb   .resize(nthrds);
  velER    .resize(nthrds);
  momD     .resize(nthrds);
  keER     .resize(nthrds);
  keIR     .resize(nthrds);
  elecOvr  .resize(nthrds, 0);
  elecAcc  .resize(nthrds, 0);
  elecTot  .resize(nthrds, 0);
  collCount.resize(nthrds);

  for (auto &v : velER) v.set_capacity(bufCap);
  for (auto &v : momD ) v.set_capacity(bufCap);
  for (auto &v : keER ) v.set_capacity(bufCap);
  for (auto &v : keIR ) v.set_capacity(bufCap);

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

  labels[elec_elec  ] = "el collisions ";
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

    if (minCollFrac > 0.0) {

				// Mean fraction in trace species
				// 
      meanF[id].clear();

      // Sample cell
      //
      pCell *samp = cell->sample;

      if (samp) {

	for (auto c : samp->children) {
	  for (auto b : c.second->bods) {
				// Particle mass accumulation
	    Particle *p = c.second->Body(b);
				// Mass-weighted trace fraction
	    speciesKey k = KeyConvert(p->iattrib[use_key]).getKey();
				// Add to species bucket
	    if (meanF[id].find(k) == meanF[id].end()) meanF[id][k] = 0.0;
	    meanF[id][k] += p->mass;
	  }
	}
      } else {

	for (auto b : cell->bods) {
				// Particle mass accumulation
	  Particle *p = cell->Body(b);
				// Mass-weighted trace fraction
	  speciesKey k = KeyConvert(p->iattrib[use_key]).getKey();
				// Add to species bucket
	  if (meanF[id].find(k) == meanF[id].end()) meanF[id][k] = 0.0;
	  meanF[id][k] += p->mass;
	}
      }
				// Normalize mass-weighted fraction
				//
      std::map<unsigned short, double> spTotl;
      for (auto v : meanF[id]) {
	unsigned short Z = v.first.first;
	if (spTotl.find(Z) == spTotl.end())
	  spTotl[Z] = v.second;
	else
	  spTotl[Z] += v.second;
      }

      for (auto &s : meanF[id]) {
	s.second /= spTotl[s.first.first];
      }
    }


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
	    Cross1 = M_PI*b*b * eVel2 * ne2 * logL;
	  }
	}

	if (i1.second>1) {
	  if (i2.second==1)
	    Cross2 = elastic(i2.first, EeV*m2/dof1) * eVel1 * ne1;
	  else {
	    double b = 0.5*esu*esu*(i2.second - 1) /
	      std::max<double>(Eerg*m2/dof1, FloorEv*eV) * 1.0e7; // nm
	    b = std::min<double>(b, ips);
	    Cross2 = M_PI*b*b * eVel1 * ne1 * logL;
	  }
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
	Cross += M_PI*b*b * eVel * meanE[id] * logL;
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
	    Cross1 = M_PI*b*b * eVel2 * ne2 * logL;
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
	      Cross2 = M_PI*b*b * eVel1 * ne1 * logL;
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
	Cross += M_PI*b*b * eVel * meanE[id] * logL;
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
  double m1   = atomic_weights[Z1] * amu;
  double m2   = atomic_weights[Z2] * amu;
  double me   = atomic_weights[0 ] * amu;
  double mu0  = m1 * m2 / (m1 + m2);
  double mu1  = m1;
  double mu2  = m2;
  double vel  = cr * UserTreeDSMC::Vunit;

  double dof1   = 1.0 + ne1;
  double dof2   = 1.0 + ne2;
  
  if (NO_DOF) dof1 = dof2 = 1.0;

  // Electron velocity equipartition factors
  //
  double eVel0 = sqrt(mu0/me);
  double eVel1 = sqrt(m1/me/dof1);
  double eVel2 = sqrt(m2/me/dof2);

  if (NO_VEL) {
    eVel0 = eVel1 = eVel2 = 1.0;
  } else if (use_elec) {
    eVel0 = eVel1 = eVel2 = 0.0;
    for (unsigned i=0; i<3; i++) {
      double rvel0 = p1->dattrib[use_elec+i] - p2->dattrib[use_elec+i];
      double rvel1 = p1->dattrib[use_elec+i] - p2->vel[i];
      double rvel2 = p2->dattrib[use_elec+i] - p1->vel[i];
      eVel0 += rvel0*rvel0;
      eVel1 += rvel1*rvel1;
      eVel2 += rvel2*rvel2;
    }
    eVel0 = sqrt(eVel0) * UserTreeDSMC::Vunit;
    eVel1 = sqrt(eVel1) * UserTreeDSMC::Vunit;
    eVel2 = sqrt(eVel2) * UserTreeDSMC::Vunit;
  }

  // Available COM energy
  //
  kEi[id] = 0.5 * mu0 * vel*vel;

  if (use_elec) {
    kEe1[id] = 0.5  * me * eVel2*eVel2;
    kEe2[id] = 0.5  * me * eVel1*eVel1;
    kEee[id] = 0.25 * me * eVel0*eVel0;
  } else {
    kEe1[id] = 0.5 * mu1 * vel*vel/dof2;
    kEe2[id] = 0.5 * mu2 * vel*vel/dof1;
  }

  // These are now ratios
  //
  eVel0 /= vel;
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
  double cross00 = 0.0;

				//-------------------------------
				// Both particles neutral
				//-------------------------------
  if (C1==1 and C2==1) {
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
      cross12 = elastic(Z1, kEe2[id]) * eVel2 * ne2 * crossfac;
      dCross[id].push_back(cross12);
      dInter[id].push_back(neut_elec_1);
    }  else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C1-1) /
	std::max<double>(kEe2[id]*eV, FloorEv*eV) * 1.0e7; // nm
      b = std::min<double>(b, ips);
      cross12 = M_PI*b*b * eVel2 * ne2 * crossfac * logL;
      dCross[id].push_back(cross12);
      dInter[id].push_back(ion_elec_1);
    }
  }
    
				//-------------------------------
				// Electrons in first particle
				//-------------------------------
  if (ne1 > 0) {
    if (C2==1) {		// Neutral atom-electron scattering
      cross21 = elastic(Z2, kEe1[id]) * eVel2 * ne2 * crossfac;
      dCross[id].push_back(cross21);
      dInter[id].push_back(neut_elec_2);
    } else {			// Rutherford scattering
      double b = 0.5*esu*esu*(C2-1) /
	std::max<double>(kEe1[id]*eV, FloorEv*eV) * 1.0e7; // nm
      b = std::min<double>(b, ips);
      cross21 = M_PI*b*b * eVel2 * ne2 * crossfac * logL;
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

    CE1[id] = ch.IonList[Q1]->collExciteCross(kEe2[id], id);

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

    double DI1 = ch.IonList[Q1]->directIonCross(kEe2[id], id);
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

    std::vector<double> RE1 = ch.IonList[Q1]->radRecombCross(kEe2[id], id);
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
    double ff2 = ch.IonList[Q2]->freeFreeCross(kEe1[id], id);
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

    CE2[id] = ch.IonList[Q2]->collExciteCross(kEe1[id], id);
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
    double DI2 = ch.IonList[Q2]->directIonCross(kEe1[id], id);
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
    std::vector<double> RE2 = ch.IonList[Q2]->radRecombCross(kEe1[id], id);
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
  return (cross00 + cross12 + cross21 + sum12 + sum21) * 1e-14 / 
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
  double mu1 = m1 * me / (m1 + me);
  double mu2 = me * m2 / (me + m2);

  // Electron velocity equipartition factors
  //
  double eVel0 = sqrt(mu0/me);
  double eVel1 = sqrt(m1/me/dof1);
  double eVel2 = sqrt(m2/me/dof2);

  if (NO_VEL) {
    eVel0 = eVel1 = eVel2 = 1.0;
  } else if (use_elec) {
    eVel0 = eVel1 = eVel2 = 0.0;
    for (unsigned i=0; i<3; i++) {
      double rvel0 = p1->dattrib[use_elec+i] - p2->dattrib[use_elec+i];
      double rvel1 = p1->dattrib[use_elec+i] - p2->vel[i];
      double rvel2 = p2->dattrib[use_elec+i] - p1->vel[i];
      eVel0 += rvel0*rvel0;
      eVel1 += rvel1*rvel1;
      eVel2 += rvel2*rvel2;
    }
    eVel0 = sqrt(eVel0) * UserTreeDSMC::Vunit;
    eVel1 = sqrt(eVel1) * UserTreeDSMC::Vunit;
    eVel2 = sqrt(eVel2) * UserTreeDSMC::Vunit;

    eVel0   /= vel;		// These are now ratios
    eVel1   /= vel;
    eVel2   /= vel;
  }
    

  // Available COM energy
  //
  kEi [id] = 0.5 * mu0 * vel*vel;
  kEe1[id] = 0.5 * mu1 * vel*vel * eVel2*eVel2/dof2;
  kEe2[id] = 0.5 * mu2 * vel*vel * eVel1*eVel1/dof1;
  kEee[id] = 0.25 * me * vel*vel * eVel0*eVel0;

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

  double cross00 = 0.0;
  double cross12 = 0.0;
  double cross21 = 0.0;

				//-------------------------------
				// Both particles neutral
				//-------------------------------
  if (C1==1 and C2==1) {
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
      cross12 = M_PI*b*b * eVel2 * ne2 * crossfac * logL;
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
      cross21 = M_PI*b*b * eVel1 * ne1 * crossfac * logL;
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
  return (cross00 + cross12 + cross21 + sum12 + sum21) * 1e-14 / 
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
      crossS += M_PI*b*b * eVel * meanE[id] * crossfac * logL;

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
  
  if (SAME_INTERACT and Z1 != Z2) return 0;

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

  double kE  = 0.5*Mu*(*cr)*(*cr);

  // For tracking energy conservation (system units)
  //
  double delE  = 0.0;

  // Now that the interactions have been calculated, create the
  // normalized cross section list to pick the interaction
  //
  std::vector<double> TotalCross;
  double tCross = 0.0;
  for (size_t i = 0; i < dCross[id].size(); i++) {
    // Sanity check (mostly for debugging, NaN should never occur)
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
    }
  }

  //----------------------------
  // Which particle interacted?
  //----------------------------
  //
  // Will be 1 or 2, dependending on which ion or neutral is
  // selected for inelastic interaction.  Will be 0 if no inealistic
  // interaction is selected.
  //
  int partflag = 0;


  // Sanity check: total cross section should be positive!
  //
  if (tCross != 0) {
    // Cumulative cross-section distribution for interaction selection
    //
    std::vector<double> CDF;
    for (size_t i = 0; i < TotalCross.size(); i++) {
      // Sanity check (mostly for debugging, NaN should never occur)
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
    //
    if (DEBUG_CR and (!DEBUG_NQ or Z1 != Z2) ) {
      //
      // Output on collisions for now . . . 
      //
      if (interFlag % 100 == 4) {
	std::cout << std::setw( 8) << "index"
		  << std::setw( 8) << "flag"
		  << std::setw(14) << "cross"
		  << std::setw(14) << "prob"
		  << std::setw(14) << "cumul"
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
		    << std::setw(14) << dCross[id][i]/tCross
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

    if (interFlag == neut_neut_1) {
      std::get<0>(ctd1->nn[id])++; 
      std::get<1>(ctd1->nn[id]) += NN;
    }

    if (interFlag == neut_elec_1) {
      std::get<0>(ctd1->ne[id])++; 
      std::get<1>(ctd1->ne[id]) += NN;
    }

    if (interFlag == ion_elec_1) {
      std::get<0>(ctd1->ie[id])++; 
      std::get<1>(ctd1->ie[id]) += NN;
    }

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

      // if (use_elec<0) delE = kEe1[id];
      delE = kEe1[id];
      if (RECOMB_IP) delE += ch.IonList[lQ(Z2, C2-1)]->ip;

      p1->iattrib[use_key] = k1.updateC(--C1);
      partflag      = 1;
      std::get<0>(ctd1->RR[id])++; 
      std::get<1>(ctd1->RR[id]) += NN;
      std::get<2>(ctd1->RR[id]) += delE * NN;
    }
    
    //-------------------------
    // Particle 2 interactions
    //-------------------------

    if (interFlag == neut_neut_2) {
      std::get<0>(ctd2->nn[id])++; 
      std::get<1>(ctd2->nn[id]) += NN;
    }

    if (interFlag == neut_elec_2) {
      std::get<0>(ctd2->ne[id])++; 
      std::get<1>(ctd2->ne[id]) += NN;
    }

    if (interFlag == ion_elec_2) {
      std::get<0>(ctd2->ie[id])++; 
      std::get<1>(ctd2->ie[id]) += NN;
    }

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

      // if (use_elec<0) delE = kEe2[id];
      delE = kEe2[id];
      if (RECOMB_IP) delE += ch.IonList[lQ(Z2, C2-1)]->ip;

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
  if (NOCOOL) delE = 0.0;

  // Elastic event
  //
  if (delE<=0.0) return ret;

  // Convert energy loss to system units
  //
  delE = delE/UserTreeDSMC::Eunit;
  
  // Assign interaction energy variables
  //
  double totE=0.0, kEe=0.0;

  // -----------------
  // ENERGY DIAGNOSTIC
  // -----------------
  // Electrons from Particle 2 have interacted with atom/ion in Particle 1
  //
  if (partflag==1) {
    totE = kE;			// KE + internal
    if (use_Eint>=0) totE += p2->dattrib[use_Eint];

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

  // Mass per particle in amu for this interaction
  //
  double m1 = atomic_weights[Z1];
  double m2 = atomic_weights[Z2];

  // Assign electron mass to doner ion particle and compute relative
  // velocity
  //
  std::vector<double> vrel(3), vcom(3), v1(3), v2(3), vcomE(3);
  double vi2 = 0.0, vf2 = 0.0;

  if (use_elec and interFlag > 100 and interFlag < 200) {

    m2 = atomic_weights[0];	// Particle 2 is the electron

    for (int k=0; k<3; k++) {
      v1[k] = p1->vel[k];	// Particle 1 is the ion
      v2[k] = p2->dattrib[use_elec+k];
      vi2  += v2[k] * v2[k];
    }

    // Secondary electron-ion scattering
    //
    if (SECONDARY_SCATTER) {
      double M1 = atomic_weights[Z2];
      double M2 = atomic_weights[ 0];
      double Mt = M1 + M2;

      for (int k=0; k<3; k++)
	vcomE[k] = (M1*p2->vel[k] + M2*p2->dattrib[use_elec+k])/Mt;
    }

  } else if (use_elec and interFlag > 200 and interFlag < 300) {

    m1 = atomic_weights[0];	// Particle 1 is the electron

    for (int k=0; k<3; k++) {
      v1[k] = p1->dattrib[use_elec+k];
      v2[k] = p2->vel[k];	// Particle 2 is the ion
      vi2  += v1[k] * v1[k];
    }

    // Secondary electron-ion scattering
    //
    if (SECONDARY_SCATTER) {
      double M1 = atomic_weights[Z1];
      double M2 = atomic_weights[ 0];
      double Mt = M1 + M2;

      for (int k=0; k<3; k++)
	vcomE[k] = (M1*p1->vel[k] + M2*p1->dattrib[use_elec+k])/Mt;
    }

  } else {
				// Neutrals or ions and electrons
    for (int k=0; k<3; k++) {
      v1[k] = p1->vel[k];
      v2[k] = p2->vel[k];
    }
  }

  // Available center of mass energy in the ballistic collision
  // (system units)
  //
  kE = 0.0;
  for (unsigned k=0; k<3; k++) {
    vcom[k] = (m1*v1[k] + m2*v2[k]) / Mt;
    kE += (v1[k] - v2[k])*(v1[k] - v2[k]);
  }

  // Relative velocity, system units
  //
  double vi = sqrt(kE);

  // Available KE in COM frame, system units
  //
  kE *= 0.5*NN*Mu;

  // Warn if energy lost is greater than total energy available to
  // lose
  //
  if (frost_warning && delE > totE)
      std::cout << "delE > KE!! (" << delE << " > " << totE
		<< "), Interaction type = " << interFlag 
		<< " kEe  = "  << kEe
		<< std::endl;

  
  // Cooling rate diagnostic histogram
  //
  if (TSDIAG && delE>0.0) {
				// Histogram index
    int indx = (int)floor(log(totE/delE)/(log(2.0)*TSPOW) + 5);
				// Floor and ceiling
    if (indx<0 ) indx = 0;
    if (indx>10) indx = 10;
				// Add entry
    EoverT[id][indx] += Mt;
  }
  
  // Time step "cooling" diagnostic
  //
  if (use_delt>=0 && delE>0.0 && totE>0.0) {
    double dtE = totE/delE * spTau[id];
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
  
  // Sufficient energy available for selected loss
  //
  if (totE > delE) {

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
				// All available energy will be lost
    lostSoFar[id] += totE;
    decolT[id]    += totE - delE;

    (*cr)          = 0.0;
    ret            = 1;		// Set error flag
    
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
  
  double cos_th = 1.0 - 2.0*(*unit)();       // Cosine and sine of
  double sin_th = sqrt(1.0 - cos_th*cos_th); // Collision angle theta
  double phi    = 2.0*M_PI*(*unit)();	     // Collision angle phi
  
  vrel[0] = vi * cos_th;	  // Compute post-collision relative
  vrel[1] = vi * sin_th*cos(phi); // velocity for an elastic 
  vrel[2] = vi * sin_th*sin(phi); // interaction
  //        ^
  //        |
  //        +---- velocity in center of mass, computed from v1, v2
  //

  // Compute the change of energy in the collision frame by computing
  // the velocity reduction factor
  //
  double vfac = 1.0;
  if (kE>0.0) vfac = totE>0.0 ? sqrt(totE/kE) : 0.0;


  // Update post-collision velocities.  In the electron version, the
  // momentum is assumed to be coupled to the ions, so the ion
  // momentum must be conserved.
  // 
  for (size_t k=0; k<3; k++) {
    v1[k] = vcom[k] + m2/Mt*vrel[k]*vfac;
    v2[k] = vcom[k] - m1/Mt*vrel[k]*vfac;
  }

  // Update electron velocties.  Electron velocity is computed so that
  // momentum is conserved ignoring the doner ion.  Use of reduction
  // factor keeps electrons and ions in equipartition.
  //
  if (use_elec and interFlag > 100 and interFlag < 200) {

    if (equiptn) {
      for (size_t k=0; k<3; k++) {
	vcom[k] = (m1*p1->vel[k] + m2*p2->dattrib[use_elec+k])/Mt;
	vrel[k] = (vcom[k] - v2[k])*Mt/m1;
      }

      for (size_t k=0; k<3; k++) {
	p1->vel[k] = vcom[k] + m2/Mt * vrel[k];
      }
    }

    // Electron from particle #2
    //
    for (size_t k=0; k<3; k++) {
      p1->vel[k] = v1[k];
      p2->dattrib[use_elec+k] = v2[k];
      vf2 += v2[k] * v2[k];
    }

    // Debug electron energy loss/gain
    //
    velER[id].push_back(vf2/vi2);
  

    // Secondary electron-ion scattering
    //
    if (SECONDARY_SCATTER) {
      double M1 = atomic_weights[Z2];
      double M2 = atomic_weights[ 0];

      for (int k=0; k<3; k++)
	  p2->vel[k] = vcomE[k] + M2/M1*(vcomE[k] - v2[k]);
    }

  } else if (use_elec and interFlag > 200 and interFlag < 300) {

    if (equiptn) {
      for (size_t k=0; k<3; k++) {
	vcom[k] = (m1*p1->dattrib[use_elec+k] + m2*p2->vel[k])/Mt;
	vrel[k] = (vcom[k] - v1[k])*Mt/m2;
      }

      for (size_t k=0; k<3; k++) {
	p2->vel[k] = vcom[k] + m1/Mt * vrel[k];
      }
    }

    // Electron from particle #1
    //
    for (size_t k=0; k<3; k++) {
      p1->dattrib[use_elec+k] = v1[k];
      p2->vel[k] = v2[k];
      vf2 += v1[k] * v1[k];
    }
    
    // Debug electron energy loss/gain
    //
    velER[id].push_back(vf2/vi2);


    // Secondary electron-ion scattering
    //
    if (SECONDARY_SCATTER) {
      double M1 = atomic_weights[Z1];
      double M2 = atomic_weights[ 0];

      for (int k=0; k<3; k++)
	p1->vel[k] = vcomE[k] + M2/M1*(vcomE[k] - v1[k]);
    }

  } else {
    for (size_t k=0; k<3; k++) {
      p1->vel[k] = v1[k];
      p2->vel[k] = v2[k];
    }
  } 

  *cr = 0.0;
  for (size_t k=0; k<3; k++) {
    double v1 = p1->vel[k];
    double v2 = p2->vel[k];
    *cr += (v1 - v2)*(v1 - v2);
  }
  *cr = sqrt(*cr);
  
  if (equiptn and use_elec) {

    if (interFlag > 100 and interFlag < 200) {

      m1 = atomic_weights[Z2];
      m2 = atomic_weights[0 ];
      Mt = m1 + m2;
      Mu = m1 * m2 / Mt;

      double KE1i = 0.0, KE2i = 0.0;
      double KE1f = 0.0, KE2f = 0.0;
      double cost = 0.0, VC2  = 0.0, VR2 = 0.0;

      for (size_t k=0; k<3; k++) {
	KE1i += p2->vel[k] * p2->vel[k];
	KE2i += p2->dattrib[use_elec+k] * p2->dattrib[use_elec+k];
	cost += p2->vel[k] * p2->dattrib[use_elec+k];

	vcom[k] = (m1*p2->vel[k] + m2*p2->dattrib[use_elec+k])/Mt;
	vrel[k] = p2->vel[k] - p2->dattrib[use_elec+k];

	VC2    += vcom[k] * vcom[k];
	VR2    += vrel[k] * vrel[k];
      }

      if (KE1i > 0.0 and KE2i > 0.0) cost /= sqrt(KE1i * KE2i);

      double dmr   = cost / (m1 - m2);
      double gamma = 1.0 + 4.0*Mt*Mu*dmr*dmr;
      double E0    = 0.5*Mt*VC2 + 0.5*Mu*VR2;

      double gamP  = 1.0 + sqrt(1.0 - 1.0/gamma);
      double gamN  = 1.0 - sqrt(1.0 - 1.0/gamma);

      double virP  = 
	(VC2 - E0/Mt*gamN)*(VC2 - E0/Mt*gamN) +
	(VR2 - E0/Mu*gamP)*(VR2 - E0/Mu*gamP) ;
	
      double virN  = 
	(VC2 - E0/Mt*gamP)*(VC2 - E0/Mt*gamP) +
	(VR2 - E0/Mu*gamN)*(VR2 - E0/Mu*gamN) ;
	
      double vcfac = 0.0, vrfac = 0.0;

      if (virP > virN) {
	vcfac = sqrt(E0/Mt*gamN);
	vrfac = sqrt(E0/Mu*gamP);
      } else {
	vcfac = sqrt(E0/Mt*gamP);
	vrfac = sqrt(E0/Mu*gamN);
      }

      if (VC2>0.0) {
	for (size_t k=0; k<3; k++) vcom[k] /= sqrt(VC2);
      } else {
	double cos_th = 1.0 - 2.0*(*unit)();
	double sin_th = sqrt(1.0 - cos_th*cos_th);
	double phi    = 2.0*M_PI*(*unit)();
	vcom[0] = cos_th;
	vcom[1] = sin_th*cos(phi);
	vcom[2] = sin_th*sin(phi);
      }

      if (VR2>0.0) {
	for (size_t k=0; k<3; k++) vrel[k] /= sqrt(VR2);
      } else {
	double cos_th = 1.0 - 2.0*(*unit)();
	double sin_th = sqrt(1.0 - cos_th*cos_th);
	double phi    = 2.0*M_PI*(*unit)();
	vrel[0] = cos_th;
	vrel[1] = sin_th*cos(phi);
	vrel[2] = sin_th*sin(phi);
      }

      for (size_t k=0; k<3; k++) {
	p2->vel[k]              = vcom[k]*vcfac + m2/Mt * vrel[k]*vrfac;
	p2->dattrib[use_elec+k] = vcom[k]*vcfac - m1/Mt * vrel[k]*vrfac;

	KE1f += p2->vel[k] * p2->vel[k];
	KE2f += p2->dattrib[use_elec+k] * p2->dattrib[use_elec+k];
      }

      KE1i *= 0.5*m1;
      KE1f *= 0.5*m1;

      KE2i *= 0.5*m2;
      KE2f *= 0.5*m2;

      double KEi = KE1i + KE2i;
      double KEf = KE1f + KE2f;

      if ( fabs(KEi - KEf) > 1.0e-14*KEi ) {
	std::cout << "Test(1): keI=[" 
		  << std::setw(16) << KE1i << ", " 
		  << std::setw(16) << KE2i << "] keF=[" 
		  << std::setw(16) << KE1f << ", " 
		  << std::setw(16) << KE2f << "] vir=[" 
		  << std::setw(16) << virP << ", "
		  << std::setw(16) << virN << "] "
		  << std::endl;
      }
    }

    if (interFlag > 200 and interFlag < 300) {

      m1 = atomic_weights[Z1];
      m2 = atomic_weights[0 ];
      Mt = m1 + m2;
      Mu = m1 * m2 / Mt;

      double KE1i = 0.0, KE2i = 0.0;
      double KE1f = 0.0, KE2f = 0.0;
      double cost = 0.0, VC2 = 0.0, VR2 = 0.0;

      for (size_t k=0; k<3; k++) {
	KE1i += p1->vel[k] * p1->vel[k];
	KE2i += p1->dattrib[use_elec+k] * p1->dattrib[use_elec+k];
	cost += p1->vel[k] * p1->dattrib[use_elec+k];

	vcom[k] = (m1*p1->vel[k] + m2*p1->dattrib[use_elec+k])/Mt;
	vrel[k] = p1->vel[k] - p1->dattrib[use_elec+k];

	VC2    += vcom[k] * vcom[k];
	VR2    += vrel[k] * vrel[k];
      }

      if (KE1i > 0.0 and KE2i > 0.0) cost /= sqrt(KE1i * KE2i);

      double dmr   = cost / (m1 - m2);
      double gamma = 1.0 + 4.0*Mt*Mu*dmr*dmr;
      double E0    = 0.5*Mt*VC2 + 0.5*Mu*VR2;

      double gamP  = 1.0 + sqrt(1.0 - 1.0/gamma);
      double gamN  = 1.0 - sqrt(1.0 - 1.0/gamma);

      double virP  = 
	(VC2 - E0/Mt*gamN)*(VC2 - E0/Mt*gamN) +
	(VR2 - E0/Mu*gamP)*(VR2 - E0/Mu*gamP) ;
	
      double virN  = 
	(VC2 - E0/Mt*gamP)*(VC2 - E0/Mt*gamP) +
	(VR2 - E0/Mu*gamN)*(VR2 - E0/Mu*gamN) ;
	
      double vcfac = 0.0, vrfac = 0.0;

      if (virP > virN) {
	vcfac = sqrt(E0/Mt*gamN);
	vrfac = sqrt(E0/Mu*gamP);
      } else {
	vcfac = sqrt(E0/Mt*gamP);
	vrfac = sqrt(E0/Mu*gamN);
      }

      if (VC2>0.0) {
	for (size_t k=0; k<3; k++) vcom[k] /= sqrt(VC2);
      } else {
	double cos_th = 1.0 - 2.0*(*unit)();
	double sin_th = sqrt(1.0 - cos_th*cos_th);
	double phi    = 2.0*M_PI*(*unit)();
	vcom[0] = cos_th;
	vcom[1] = sin_th*cos(phi);
	vcom[2] = sin_th*sin(phi);
      }

      if (VR2>0.0) {
	for (size_t k=0; k<3; k++) vrel[k] /= sqrt(VR2);
      } else {
	double cos_th = 1.0 - 2.0*(*unit)();
	double sin_th = sqrt(1.0 - cos_th*cos_th);
	double phi    = 2.0*M_PI*(*unit)();
	vrel[0] = cos_th;
	vrel[1] = sin_th*cos(phi);
	vrel[2] = sin_th*sin(phi);
      }

      for (size_t k=0; k<3; k++) {
	p1->vel[k]              = vcom[k]*vcfac + m2/Mt * vrel[k]*vrfac;
	p1->dattrib[use_elec+k] = vcom[k]*vcfac - m1/Mt * vrel[k]*vrfac;

	KE1f += p1->vel[k] * p1->vel[k];
	KE2f += p1->dattrib[use_elec+k] * p1->dattrib[use_elec+k];
      }

      KE1i *= 0.5*m1;
      KE1f *= 0.5*m1;

      KE2i *= 0.5*m2;
      KE2f *= 0.5*m2;

      double KEi = KE1i + KE2i;
      double KEf = KE1f + KE2f;

      if ( fabs(KEi - KEf) > 1.0e-14*KEi ) {
	std::cout << "Test(1): keI=[" 
		  << std::setw(16) << KE1i << ", " 
		  << std::setw(16) << KE2i << "] keF=[" 
		  << std::setw(16) << KE1f << ", " 
		  << std::setw(16) << KE2f << "] vir=[" 
		  << std::setw(16) << virP << ", "
		  << std::setw(16) << virN << "] "
		  << std::endl;
      }
    }
  }

  // Scatter electrons
  //
  if (esType == always and C1>1 and C2>1) {
    double vi = 0.0;
    for (int k=0; k<3; k++) {
      double d1 = p1->dattrib[use_elec+k];
      double d2 = p2->dattrib[use_elec+k];
      vcom[k] = 0.5*(d1 + d2);
      vi     += (d1 - d2) * (d1 - d2);
    }
    vi = sqrt(vi);

    double cos_th = 1.0 - 2.0*(*unit)();       // Cosine and sine of
    double sin_th = sqrt(1.0 - cos_th*cos_th); // Collision angle theta
    double phi    = 2.0*M_PI*(*unit)();	     // Collision angle phi
    
    vrel[0] = vi * cos_th;	  // Compute post-collision relative
    vrel[1] = vi * sin_th*cos(phi); // velocity for an elastic 
    vrel[2] = vi * sin_th*sin(phi); // interaction

    for (int k=0; k<3; k++) {
      p1->dattrib[use_elec+k] = vcom[k] + 0.5*vrel[k];
      p2->dattrib[use_elec+k] = vcom[k] - 0.5*vrel[k];
    }

    return 0;
  }

  return ret;
}

int CollideIon::inelasticWeight(pCell* const c, 
				Particle* const _p1, Particle* const _p2,
				double *cr, int id)
{
  int ret       =  0;		// No error (flag)
  int interFlag = -1;		// Invalid value by default

  Particle* p1  = _p1;		// Copy pointers for swapping, if
  Particle* p2  = _p2;		// necessary


  // Species keys for pointers before swapping
  //
  KeyConvert k1(p1->iattrib[use_key]);
  KeyConvert k2(p2->iattrib[use_key]);

  collTDPtr ctd1 = (*collD)[k1.getKey()];
  collTDPtr ctd2 = (*collD)[k2.getKey()];

  unsigned short Z1 = k1.getKey().first, C1 = k1.getKey().second;
  unsigned short Z2 = k2.getKey().first, C2 = k2.getKey().second;

  if (SAME_INTERACT and Z1 != Z2) return 0;

  // Particle 1 is assumed to be the "dominant" species and Particle 2
  // is assumed to be the "trace" species (or another "dominant").
  // Swap particle pointers if necessary.
  //
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
      
  // Debugging test
  //
  double p1E = 0.0, p2E = 0.0;
  if (use_cons_test and use_cons>=0) {
    p1E = p1->dattrib[use_cons];
    p2E = p2->dattrib[use_cons];
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
  for (size_t i = 0; i < dCross[id].size(); i++) {

    // Sanity check (mostly for debugging, NaN should never occur)
    //
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

      bool ok = false;		// Reject all interactions by default

      // Accumulate the list here
      //
      if (NoDelC)  {
	ok = true;
				// Pass events that are NOT ionization
				// or recombination, or both
	if (NoDelC & 0x1 and dInter[id][i] % 100 == recomb) ok = false;
	if (NoDelC & 0x2 and dInter[id][i] % 100 == ionize) ok = false;

      } else if (scatter) {
				// Only pass elastic scattering events
	if (dInter[id][i] % 100 < 3) ok = true;

				// Otherwise, test all events . . . 
      } else {
				// Test for Particle #1 collisional excitation
	if (dInter[id][i] == 104) {
	  double frac = meanF[id][k1.getKey()];
	  if (frac > minCollFrac) {
	    ok = true;
	  }
	}
				// Test for Particle #2 collisional excitation
	else if (dInter[id][i] == 204) {
	  double frac = meanF[id][k2.getKey()];
	  if (frac > minCollFrac) {
	    ok = true;
	  }
	}
	else {			// Pass all other interactions . . . 
	  ok = true;
	}
      }

      if (ok) tCross += dCross[id][i];

      TotalCross.push_back(tCross);
    }
  }

  //
  // Cross section scale factor
  //
  double scaleCrossSection = tCross/csections[id][k1.getKey()][k2.getKey()] *
    1e-14 / (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);

  NN *= scaleCrossSection;

  //----------------------------
  // Which particle interacted?
  //----------------------------
  //
  // Will be 1 or 2, dependending on which ion or neutral is
  // selected for inelastic interaction.  Will be 0 if no inealistic
  // interaction is selected.
  //
  int partflag = 0;

  // NOCOOL debugging
  //
  double NCXTRA = 0.0;

  double iI1 = 0.0, iE1 = 0.0;
  double iI2 = 0.0, iE2 = 0.0;

  if (NOCOOL) {
    for (size_t k=0; k<3; k++) {
      iI1 += p1->vel[k]*p1->vel[k];
      iI2 += p2->vel[k]*p2->vel[k];
      if (use_elec>=0) {
	if (C1>1) iE1 += p1->dattrib[use_elec+k]*p1->dattrib[use_elec+k];
	if (C2>1) iE2 += p2->dattrib[use_elec+k]*p2->dattrib[use_elec+k];
      }
    }
    iI1 *= 0.5*p1->mass;
    iI2 *= 0.5*p2->mass;
    iE1 *= 0.5*p1->mass*atomic_weights[0]/atomic_weights[Z1];
    iE2 *= 0.5*p2->mass*atomic_weights[0]/atomic_weights[Z2];
  }

  // Sanity check: total cross section should be positive!
  //
  if (tCross > 0.0) {

    // Cumulative cross-section distribution for interaction selection
    //
    std::vector<double> CDF;
    for (size_t i = 0; i < TotalCross.size(); i++) {
      // Sanity check (mostly for debugging, NaN should never occur)
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
    //
    if (DEBUG_CR and (!DEBUG_NQ or Z1 != Z2) ) {
      speciesKey i1 = k1.getKey();
      speciesKey i2 = k2.getKey();
      double cfac = 1e-14 / (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
      //
      // Output on collisions for now . . . 
      //
      std::cout << std::setw( 8) << "index"
		<< std::setw( 8) << "flag"
		<< std::setw(14) << "cross"
		<< std::setw(14) << "prob"
		<< std::setw(14) << "cumul"
		<< std::setw(14) << "tCross/max"
		<< std::setw(18) << "type label"
		<< std::endl
		<< std::setw( 8) << "-----"
		<< std::setw( 8) << "-----"
		<< std::setw(14) << "---------"
		<< std::setw(14) << "---------"
		<< std::setw(14) << "---------"
		<< std::setw(14) << "---------"
		<< std::setw(18) << "---------------"
		<< std::endl;
      for (size_t i = 0; i < dCross[id].size(); i++) {
	std::cout << std::setw( 8) << i
		  << std::setw( 8) << dInter[id][i]
		  << std::setw(14) << dCross[id][i]
		  << std::setw(14) << dCross[id][i]/tCross
		  << std::setw(14) << CDF[i]
		  << std::setw(14) << dCross[id][i]/csections[id][i1][i2] * cfac
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

    if (interFlag == neut_neut_1) {
      std::get<0>(ctd1->nn[id])++; 
      std::get<1>(ctd1->nn[id]) += NN;
    }

    if (interFlag == neut_elec_1) {
      std::get<0>(ctd1->ne[id])++; 
      std::get<1>(ctd1->ne[id]) += NN;
    }

    if (interFlag == ion_elec_1) {
      std::get<0>(ctd1->ie[id])++; 
      std::get<1>(ctd1->ie[id]) += NN;
    }

    if (interFlag == free_free_1) {
      delE          = IS.selectFFInteract(ch.IonList[Q1], id);
      partflag      = 1;
      if (NO_FF_E) delE = 0.0;
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
      if (NO_ION_E) delE = 0.0;
      std::get<0>(ctd1->CI[id])++; 
      std::get<1>(ctd1->CI[id]) += Wb;
      std::get<2>(ctd1->CI[id]) += delE * NN;
    }

    if (interFlag == recomb_1) {

      p1->iattrib[use_key] = k1.updateC(--C1);
      partflag      = 1;

      // if (use_elec<0) delE = kEe1[id];
      delE = kEe1[id];
      if (RECOMB_IP) delE += ch.IonList[lQ(Z1, C1)]->ip;

      std::get<0>(ctd1->RR[id])++; 
      std::get<1>(ctd1->RR[id]) += Wb;
      std::get<2>(ctd1->RR[id]) += delE * NN;

      // Add the KE from the recombined electron back to the free pool
      //
      if (NOCOOL and !NOCOOL_ELEC and C1==1 and use_cons>=0) {
	double lKE = 0.0, fE = 0.5*Wa*atomic_weights[0];
	for (size_t k=0; k<3; k++) {
	  double t = p1->dattrib[use_elec+k];
	  lKE += fE*t*t;
	}

	NCXTRA += lKE;

	if (q<1)
	  p1->dattrib[use_cons] += lKE;
	else {
	  p1->dattrib[use_cons] += lKE * 0.5;
	  p2->dattrib[use_cons] += lKE * 0.5;
	}
      }
    }
    
    //-------------------------
    // Particle 2 interactions
    //-------------------------

    if (interFlag == neut_neut_2) {
      std::get<0>(ctd2->nn[id])++; 
      std::get<1>(ctd2->nn[id]) += NN;
    }

    if (interFlag == neut_elec_2) {
      std::get<0>(ctd2->ne[id])++; 
      std::get<1>(ctd2->ne[id]) += NN;
    }

    if (interFlag == ion_elec_2) {
      std::get<0>(ctd2->ie[id])++; 
      std::get<1>(ctd2->ie[id]) += NN;
    }

    if (interFlag == free_free_2) {
      delE          = IS.selectFFInteract(ch.IonList[Q2], id);
      partflag      = 2;
      if (NO_FF_E) delE = 0.0;
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
      partflag     = 2;
      if (NO_ION_E) delE = 0.0;
      std::get<0>(ctd2->CI[id])++; 
      std::get<1>(ctd2->CI[id]) += Wb;
      std::get<2>(ctd2->CI[id]) += delE * NN;
    }

    if (interFlag == recomb_2) {

      p2->iattrib[use_key] = k2.updateC(--C2);
      partflag     = 2;

      // if (use_elec<0) delE = kEe2[id];
      delE = kEe2[id];
      if (RECOMB_IP) delE += ch.IonList[lQ(Z2, C2)]->ip;

      std::get<0>(ctd2->RR[id])++; 
      std::get<1>(ctd2->RR[id]) += Wb;
      std::get<2>(ctd2->RR[id]) += delE * NN;

      // Add the KE from the recombined electron back to the free pool
      //
      if (NOCOOL and !NOCOOL_ELEC and C2==1 and use_cons>=0) {
	double lKE = 0.0, fE = 0.5*Wb*atomic_weights[0];
	for (size_t k=0; k<3; k++) {
	  double t = p2->dattrib[use_elec+k];
	  lKE += fE*t*t;
	}

	NCXTRA += lKE;

	if (q<1)
	  p1->dattrib[use_cons] += lKE;
	else {
	  p1->dattrib[use_cons] += lKE * 0.5;
	  p2->dattrib[use_cons] += lKE * 0.5;
	}
      }
    }

    // Convert to super particle
    //
    delE *= NN;
    
    // Convert back to cgs
    //
    delE = delE * eV;
  }
  
  // Collision counts
  //
  if (COLL_SPECIES) {
    dKey dk(k1.getKey(), k2.getKey());
    if (collCount[id].find(dk) == collCount[id].end()) collCount[id][dk] = ccZ;
    if (interFlag % 100 <= 2) collCount[id][dk][0]++;
    else                      collCount[id][dk][1]++;
  }

  // Debugging test
  //
  if (SAME_IONS_SCAT and interFlag % 100 <= 2) {
    if (Z1 != Z2) return 0;
  }

  // Work vectors
  //
  std::vector<double> vrel(3), vcom(3), v1(3), v2(3);

  // For elastic interactions, delE == 0
  //
  assert(delE >= 0.0);

  // Artifically prevent cooling by setting the energy removed from
  // the COM frame to zero
  //
  if (NOCOOL) delE = 0.0;

  // Convert energy loss to system units
  //
  delE = delE/UserTreeDSMC::Eunit;

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

  //
  // Perform energy adjustment in ion, system COM frame with system
  // mass units
  //

  // Mass per particle in amu for this interaction
  //
  double m1 = atomic_weights[Z1];
  double m2 = atomic_weights[Z2];

  // Assign electron mass to doner ion particle and compute relative
  // velocity
  //
  double vi2 = 0.0, vf2 = 0.0;

  if (use_elec and interFlag > 100 and interFlag < 200) {

    m2 = atomic_weights[0];	// Particle 2 is the electron

    for (int k=0; k<3; k++) {
      v1[k] = p1->vel[k];	// Particle 1 is the ion
      v2[k] = p2->dattrib[use_elec+k];
      vi2  += v2[k] * v2[k];
    }

  } else if (use_elec and interFlag > 200 and interFlag < 300) {

    m1 = atomic_weights[0];	// Particle 1 is the electron

    for (int k=0; k<3; k++) {
      v1[k] = p1->dattrib[use_elec+k];
      v2[k] = p2->vel[k];	// Particle 2 is the ion
      vi2  += v1[k] * v1[k];
    }

  } else {
				// Neutrals or ions and electrons
    for (int k=0; k<3; k++) {
      v1[k] = p1->vel[k];
      v2[k] = p2->vel[k];
    }
  }

  // For debugging kinetic energy bookkeeping
  //
  double KE1i = 0.0, KE2i = 0.0;
  double KE1f = 0.0, KE2f = 0.0;

  if (KE_DEBUG) {
    for (auto v : v1) KE1i += v*v;
    for (auto v : v2) KE2i += v*v;
  }

  // Total effective mass in the collision (atomic mass units)
  //
  double Mt = m1 + m2;

  // Reduced mass (atomic mass units)
  //
  double Mu = m1 * m2 / Mt;


  // Available center of mass energy in the ballistic collision
  // (system units)
  //
  double kE = 0.0;
  for (unsigned k=0; k<3; k++) {
    vcom[k] = (m1*v1[k] + m2*v2[k]) / Mt;
    kE += (v1[k] - v2[k])*(v1[k] - v2[k]);
  }

  // Relative velocity, system units
  //
  double vi = sqrt(kE);

  // Available KE in COM frame, system units
  //
  kE *= 0.5*Wa*q*Mu;

  // Total energy available in COM after removing radiative and
  // collisional loss.  A negative value for totE will be handled
  // below . . .
  //
  double totE  = kE - delE;

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

  double Exs = 0.0;		// Exs variable for KE debugging only

  bool in_exactE = false;
  
  if (use_cons >= 0) {

    if (SAME_TRACE_SUPP) {
      //
      // Override special trace species treatment
      //
      double del = 0.0;

      if (use_elec>=0) {

				// Particle 1: ion
				// Particle 2: electron
	if (interFlag > 100 and interFlag < 200) {

	  del += p1->dattrib[use_cons];
	  p1->dattrib[use_cons] = 0.0;

	  del += p2->dattrib[use_elec+3];
	  p2->dattrib[use_elec+3] = 0.0;

				// Particle 1: electron
				// Particle 2: ion
	} else if (interFlag > 200 and interFlag < 300) {
	  
	  del += p1->dattrib[use_elec+3];
	  p1->dattrib[use_elec+3] = 0.0;

	  del += p2->dattrib[use_cons];
	  p2->dattrib[use_cons] = 0.0;

				// Particle 1: ion
				// Particle 2: ion
	} else {

	  del += p1->dattrib[use_cons];
	  p1->dattrib[use_cons] = 0.0;

	  del += p2->dattrib[use_cons];
	  p2->dattrib[use_cons] = 0.0;

	}

      } else {
	del = p1->dattrib[use_cons] + p2->dattrib[use_cons];
	p1->dattrib[use_cons] = p2->dattrib[use_cons] = 0.0;
      }

      Exs  += del;
      totE += del;

    } else {

      //
      // Not a trace interaction
      //
      if (Z1 == Z2) {
	
	double del = 0.0;

	if (use_elec>=0) {
				// Particle 1: ion
				// Particle 2: electron
	  if (interFlag > 100 and interFlag < 200) {

	    del += p1->dattrib[use_cons];
	    p1->dattrib[use_cons] = 0.0;

	    del += p2->dattrib[use_elec+3];
	    p2->dattrib[use_elec+3] = 0.0;

				// Particle 1: electron
				// Particle 2: ion
	  } else if (interFlag > 200 and interFlag < 300) {
	  
	    del += p1->dattrib[use_elec+3];
	    p1->dattrib[use_elec+3] = 0.0;
	    
	    del += p2->dattrib[use_cons];
	    p2->dattrib[use_cons] = 0.0;

				// Particle 1: ion
				// Particle 2: ion
	  } else {

	    del += p1->dattrib[use_cons];
	    p1->dattrib[use_cons] = 0.0;

	    del += p2->dattrib[use_cons];
	    p2->dattrib[use_cons] = 0.0;
	    
	  }

	} else {
	  del = p1->dattrib[use_cons] + p2->dattrib[use_cons];
	  p1->dattrib[use_cons] = p2->dattrib[use_cons] = 0.0;
	}


	Exs  += del;
	totE += del;
	
      } else if (ExactE) {
	
	// Particle 1: ion
	// Particle 2: electron
	if (interFlag > 100 and interFlag < 200) {
	  p1->dattrib[use_cons]   += -0.5*delE;
	  p2->dattrib[use_elec+3] += -0.5*delE;
	}
	// Particle 1: electron
	// Particle 2: ion
	else if (interFlag > 200 and interFlag < 300) {
	  p1->dattrib[use_elec+3] += -0.5*delE;
	  p2->dattrib[use_cons  ] += -0.5*delE;
	} 
	// Neutral interaction
	else {
	  p1->dattrib[use_cons  ] += -0.5*delE;
	  p2->dattrib[use_cons  ] += -0.5*delE;
	}

	// Reset total energy to initial energy, deferring an changes to
	// non-trace interactions
	//
	totE = kE;

	in_exactE = true;
      }

    } // end : SAME_TRAC_SUPP if/then

  } // end: trace-particle energy loss assignment


  double cos_th = 1.0 - 2.0*(*unit)();       // Cosine and sine of
  double sin_th = sqrt(1.0 - cos_th*cos_th); // Collision angle theta
  double phi    = 2.0*M_PI*(*unit)();	     // Collision angle phi
  
  vrel[0] = vi * cos_th;	  // Compute post-collision relative
  vrel[1] = vi * sin_th*cos(phi); // velocity for an elastic 
  vrel[2] = vi * sin_th*sin(phi); // interaction
  //        ^
  //        |
  //        +---- velocity in center of mass, computed from v1, v2
  //

  // Attempt to defer negative energy adjustment
  //
  double missE = std::min<double>(0.0, totE);


  // Compute the change of energy in the collision frame by computing
  // the velocity reduction factor
  //
  double vfac = 1.0;
  if (kE>0.0) vfac = totE>0.0 ? sqrt(totE/kE) : 0.0;


  // Use explicit energy conservation algorithm
  //
  double vrat = 1.0;
  std::vector<double> w1(v1);

  if (ExactE and q < 1.0) {

    double v1i2 = 0.0, b1f2 = 0.0, v2i2 = 0.0, b2f2 = 0.0;
    std::vector<double> uu(3), vv(3);
    for (size_t k=0; k<3; k++) {
      uu[k] = vcom[k] + m2/Mt*vrel[k];
      vv[k] = vcom[k] - m1/Mt*vrel[k];
      v1i2 += v1[k]*v1[k];
      v2i2 += v2[k]*v2[k];
      b1f2 += uu[k]*uu[k];
      b2f2 += vv[k]*vv[k];
    }

    if (AlgOrth) {

      // Cross product to determine orthgonal direction
      //
      w1[0] = uu[1]*v1[2] - uu[2]*v1[1];
      w1[1] = uu[2]*v1[0] - uu[0]*v1[2];
      w1[2] = uu[0]*v1[1] - uu[1]*v1[0];

      // Normalize
      //
      double wnrm = 0.0;
      for (auto   v : w1) wnrm += v*v;
      for (auto & v : w1) v *= 1.0/sqrt(wnrm);
      
      vrat  = sqrt((1.0 - q)*(q*b1f2 + v1i2))/(1.0 - q);

    } else {

      double qT = 0.0;
      for (size_t k=0; k<3; k++) qT += v1[k]*uu[k];
      
      if (v1i2 > 0.0 and b1f2 > 0.0) qT *= q/sqrt(v1i2 * b1f2);
      
      double vh1f  = 
	( -sqrt(b1f2)*qT + sqrt(qT*qT*b1f2 + (1.0 - q)*(q*b1f2 + v1i2)) )/(1.0 - q);

      vrat = vh1f / sqrt(v1i2);
    }

    // Test
    double v1f2 = 0.0;
    for (size_t k=0; k<3; k++) {
      double vv = (1.0 - q)*w1[k]*vrat + q*uu[k];
      v1f2 += vv*vv;
    }
    
    double KE_1i = 0.5*m1*v1i2;
    double KE_2i = 0.5*q*m2*v2i2;
    double KE_1f = 0.5*m1*v1f2;
    double KE_2f = 0.5*q*m2*b2f2;
    
    double KEi   = KE_1i + KE_2i;
    double KEf   = KE_1f + KE_2f;
    double difE  = KEi - KEf;

    if (fabs(difE)/(KEi + KEf) > 1.0e-10) {
      std::cout << "Ooops, delE = " << difE
		<< ", totE = " << KEi + KEf << std::endl;
    }
  }

  // Update post-collision velocities.  In the electron version, the
  // momentum is assumed to be coupled to the ions, so the ion
  // momentum must be conserved.  Particle 2 is trace by construction.
  // 
  double deltaKE = 0.0, qKEfac = 0.5*Wa*m1*q*(1.0 - q);
  for (size_t k=0; k<3; k++) {
    double v0 = vcom[k] + m2/Mt*vrel[k]*vfac;

    if (!ExactE) 
      deltaKE += (v0 - v1[k])*(v0 - v1[k]) * qKEfac;

    v1[k] = (1.0 - q)*w1[k]*vrat + q*v0;
    v2[k] = vcom[k] - m1/Mt*vrel[k]*vfac;
  }

  // Save energy adjustments for next interation.  Split between like
  // species ONLY.
  //
  if (use_cons >= 0) {

    double del = missE;

    // Energy is added to electron KE for use_elec >= 0
    if (use_elec < 0) 
      del += deltaKE;

    else if (C1==1 and C2==1)
      del += deltaKE;
      
    if (SAME_TRACE_SUPP) {
      //
      // Override default trace species treatment; split energy
      // adjustment between interaction particles
      //
      if (C1 == 1 or use_elec<0)
	p1->dattrib[use_cons  ] += 0.5*del;
      else
	p1->dattrib[use_elec+3] += 0.5*del;

      if (C1 == 2 or use_elec<0)
	p2->dattrib[use_cons]   += 0.5*del;
      else
	p2->dattrib[use_elec+3] += 0.5*del;

    } else {
      //
      // Split energy adjustment between like species ONLY.
      // Otherwise, assign to non-trace particle.
      //
      if (Z1 == Z2) {		
	p1->dattrib[use_cons] += 0.5*del;
	p2->dattrib[use_cons] += 0.5*del;
      } else {
	if (C1 == 1 or use_elec<0)
	  p1->dattrib[use_cons  ] += del;
	else
	  p1->dattrib[use_elec+3] += del;
      }
    }
  }

  // Update particle velocties
  // -------------------------
  //
  // Electron velocity is computed so that momentum is conserved
  // ignoring the doner ion
  //
  bool electronic = false;

  if (use_elec and interFlag > 100 and interFlag < 200) {
    
    electronic = true;
    
    if (equiptn) {
      for (size_t k=0; k<3; k++) {
	vcom[k] = (m1*p1->vel[k] + m2*p2->dattrib[use_elec+k])/Mt;
	vrel[k] = (vcom[k] - v2[k])*Mt/m1;
      }

      for (size_t k=0; k<3; k++) {
	p1->vel[k] = vcom[k] + m2/Mt * vrel[k];
      }
    }

    // Upscale electron velocity to conserve energy
    //
    double vfac = 1.0;
    if (Z1 != Z2) {
      if (TRACE_ELEC and !TRACE_REAPPLY) {
	p1->dattrib[use_cons  ] += deltaKE * (1.0 - TRACE_FRAC);
	p2->dattrib[use_elec+3] += deltaKE * TRACE_FRAC;
      } else {
	double ke2 = 0.0;
	for (auto v : v2) ke2 += v*v;
	ke2 *= 0.5*Wb*m2;
	vfac = sqrt(1.0 + deltaKE/ke2);
      }
    }
    
    // Electron from particle #2
    //
    double elecE0 = 0.0;
    for (size_t k=0; k<3; k++) {
				// Compute energy
      v2[k] *= vfac;
      vf2 += v2[k] * v2[k];
				// Assign to particles
      p1->vel[k] = v1[k];
      elecE0 += v2[k] * v2[k];
      p2->dattrib[use_elec+k] = v2[k];
				// Zero velocity for recombined
				// electron
      if (interFlag == recomb_1 and C1==1)
	p1->dattrib[use_elec+k] = 0.0;
    }
    
				// Duplicate electron energy if previously
				// neutral
				//
    if (CLONE_ELEC and interFlag == ionize_1 and C1==2) {

      double EE1 = 0.0, EE2 = 0.0;
	for (size_t k=0; k<3; k++) {
	  EE1 += v1[k]*v1[k];
	  EE2 += v2[k]*v2[k];
      }

      for (size_t k=0; k<3; k++)
	p2->dattrib[use_elec+k] = v2[k] * EE1/EE2;

				// Share electron energy if previously
				// neutral
				//
    } else if (!NOSHARE_ELEC and interFlag == ionize_1 and C1==2) {
				// KE prefactor for Particle #1
      double fE1 = 0.5*Wa*atomic_weights[0]; 
				// KE prefactor for Particle #2
      double fE2 = 0.5*Wb*atomic_weights[0];
				// For energy conservation
      double elecE1 = 0.0, elecE2 = 0.0;

      //
      // Split the energy between the outgoing electrons
      //

      //
      // Initial electron energy of Particle #1 is 0 (i.e. neutral)
      // Total electron energy is that of Particle #2:
      // E0= 1/2*m_e*Wb*v0^2
      //
      // Split between two electrons:
      // E1 = (1 - u)*E0 = (1 - u)*1/2*m_e*Wb*v0^2 = 1/2*m_e*Wa*v1^2
      // ---> v1^2 = (1-u)*q*v0^2
      //
      // E2 = u*E0 = u*1/2*m_e*Wb*v0^2 = 1/2*m_e*Wb*v2^2
      // ---> v2^2 = u*v0^2
      //
      double u  = (*unit)();
      double vs1 = sqrt(q*(1.0 - u));
      double vs2 = sqrt(u);

      for (size_t k=0; k<3; k++) {
	double t1 = vs1*v2[k];
	double t2 = vs2*v2[k];
	elecE1 += fE1*t1*t1;
	elecE2 += fE2*t2*t2;
	p1->dattrib[use_elec+k] = t1;
	p2->dattrib[use_elec+k] = t2;
      }

      elecE0 *= fE2;

      double deltaE_e = elecE0 - elecE1 - elecE2;
      
      NCXTRA += deltaE_e;

      if (use_cons>=0) {
	if (q<1.0) {
	  p1->dattrib[use_elec+3] += deltaE_e;
	} else {
	  p1->dattrib[use_elec+3] += 0.5*deltaE_e;
	  p2->dattrib[use_elec+3] += 0.5*deltaE_e;
	}
      }
    }

    // For diagnostic electron energy loss/gain distribution
    //
    velER[id].push_back(vf2/vi2);
    
  } else if (use_elec and interFlag > 200 and interFlag < 300) {

    electronic = true;

    if (equiptn) {
      for (size_t k=0; k<3; k++) {
	vcom[k] = (m1*p1->dattrib[use_elec+k] + m2*p2->vel[k])/Mt;
	vrel[k] = (vcom[k] - v1[k])*Mt/m2;
      }

      for (size_t k=0; k<3; k++) {
	p2->vel[k] = vcom[k] + m1/Mt * vrel[k];
      }
    }

    // Upscale electron velocity to conserve energy
    //
    double vfac = 1.0;
    if (Z1 != Z2) {
      if (TRACE_ELEC and !TRACE_REAPPLY) {
	p1->dattrib[use_elec+3] += deltaKE * TRACE_FRAC;
	p2->dattrib[use_cons]   += deltaKE * (1.0 - TRACE_FRAC);
      } else {
	double ke1 = 0.0;
	for (auto v : v1) ke1 += v*v;
	ke1 *= 0.5*Wa*m1;
	vfac = sqrt(1.0 + deltaKE/ke1);
      }
    }

    // Electron from particle #1
    //
    double elecE0 = 0.0;
    for (size_t k=0; k<3; k++) {
				// Compute energy
      v1[k] *= vfac;
      vf2 += v1[k] * v1[k];
				// Assign to particles
      elecE0 += v1[k] * v1[k];
      p1->dattrib[use_elec+k] = v1[k];
      p2->vel[k] = v2[k];
				// Zero velocity for recombined
				// electron
      if (interFlag == recomb_2 and C2==1)
	p2->dattrib[use_elec+k] = 0.0;
    }
    
				// Duplicate electron energy if previously
				// neutral
				//
    if (CLONE_ELEC and interFlag == ionize_2 and C2==2) {

      double EE1 = 0.0, EE2 = 0.0;
	for (size_t k=0; k<3; k++) {
	  EE1 += v1[k]*v1[k];
	  EE2 += v2[k]*v2[k];
      }

      for (size_t k=0; k<3; k++)
	p1->dattrib[use_elec+k] = v2[k] * EE2/EE1;

				// Share electron energy if previously
				// neutral
				//
    } else if (!NOSHARE_ELEC and interFlag == ionize_2 and C2==2) {
				// KE prefactor for Particle #1
      double fE1 = 0.5*Wa*atomic_weights[0]; 
				// KE prefactor for Particle #2
      double fE2 = 0.5*Wb*atomic_weights[0];
				// For energy conservation
      double elecE1 = 0.0, elecE2 = 0.0;

      //
      // Split fraction q of the energy between the outgoing electrons
      //

      //
      // Initial electron energy of Particle #2 is 0 (neutral)
      // Total electron energy is that of Particle #1:
      // E0= 1/2*m_e*Wa*v0^2
      //
      // Split initial particle by q, then between two electrons:
      // E1 = (1 - q)*E0 + q*(1 - u)*E0
      //    = [(1 - q) + q*(1 - u)]*1/2*m_e*Wa*v0^2 = 1/2*m_e*Wa*v1^2
      // ---> v1^2 = [1-q + q*(1-u)]*v0^2
      //
      // E2 = u*q*E0 = u*1/2*m_e*Wb*v0^2 = 1/2*m_e*Wb*v2^2
      // ---> v2^2 = u*v0^2
      //

      double u  = (*unit)();
      double vs1 = sqrt((1.0 - q) + q*(1.0 - u));
      double vs2 = sqrt(u);

      for (size_t k=0; k<3; k++) {
	double t1 = vs1*v1[k];
	double t2 = vs2*v1[k];
	elecE1 += fE1*t1*t1;
	elecE2 += fE2*t2*t2;
	p1->dattrib[use_elec+k] = t1;
	p2->dattrib[use_elec+k] = t2;
      }

      elecE0 *= fE1;

      double deltaE_e = elecE0 - elecE1 - elecE2;
      
      NCXTRA += deltaE_e;

      if (use_cons>=0) {
	if (q<1.0) {
	  p1->dattrib[use_elec+3] += deltaE_e;
	} else {
	  p1->dattrib[use_elec+3] += 0.5*deltaE_e;
	  p2->dattrib[use_elec+3] += 0.5*deltaE_e;
	}
      }
    }

    // For diagnostic electron energy loss/gain distribution
    //
    velER[id].push_back(vf2/vi2);

  } else {
    for (size_t k=0; k<3; k++) {
      p1->vel[k] = v1[k];
      p2->vel[k] = v2[k];
    }
  } 

  // KE debugging
  //
  if (KE_DEBUG) {

    for (auto v : v1) KE1f += v*v;
    for (auto v : v2) KE2f += v*v;

				// Pre collision KE
    KE1i *= 0.5*Wa*m1;
    KE2i *= 0.5*Wb*m2;
				// Post collision KE
    KE1f *= 0.5*Wa*m1;
    KE2f *= 0.5*Wb*m2;

    double tKEi = KE1i + KE2i;	// Total pre collision KE
    double tKEf = KE1f + KE2f;	// Total post collision KE
    double dKE  = tKEi - tKEf;	// Energy balance

    if (m1<1.0) {
      if (KE1i > 0) keER[id].push_back((KE1i - KE1f)/KE1i);
      if (KE2i > 0) keIR[id].push_back((KE2i - KE2f)/KE2i);
    } 

    if (m2<1.0) {
      if (KE1i > 0) keIR[id].push_back((KE1i - KE1f)/KE1i);
      if (KE2i > 0) keER[id].push_back((KE2i - KE2f)/KE2i);
    }
				// Check energy balance including excess
    double testE = dKE;

    if (Z1 == Z2 or !ExactE) testE -= delE + missE;

				// Add in energy loss/gain
    if (Z1==Z2 or SAME_TRACE_SUPP)
      testE += Exs;
				// Correct for trace-algorithm excess
    else if ( (C1==1 and C2==1) or (electronic and !TRACE_REAPPLY) )
      testE -= deltaKE;

    if (fabs(testE) > DEBUG_THRESH*(tKEi+tKEf) )
      std::cout << "Total ("<< m1 << "," << m2 << ") = " 
		<< std::setw(14) << testE
		<< ", dKE=" << std::setw(14) << dKE
		<< ", KE0=" << std::setw(14) << kE
		<< ", tot=" << std::setw(14) << totE
		<< ", com=" << std::setw(14) << deltaKE
		<< ", exs=" << std::setw(14) << Exs
		<< ", del=" << std::setw(14) << delE
		<< ", mis=" << std::setw(14) << missE
		<< ", fac=" << std::setw(14) << vfac
		<< (in_exactE ? ", in ExactE" : "")
		<< std::endl;

  } // Energy conservation debugging diagnostic (KE_DEBUG)
  

  // Debugging test
  //
  if (use_cons_test and use_cons>=0) {

    if (p1->dattrib[use_cons] - p1E > 0.0 and p1->dattrib[use_cons]>0.0) {
      std::cout << "P2 above zero: dif=" 
		<< p1->dattrib[use_cons] - p1E          << std::endl
		<< "    x_f="  << p1->dattrib[use_cons] << std::endl
		<< "    x_i="  << p1E                   << std::endl
		<< "    exs="  << Exs                   << std::endl
		<< "   misE="  << missE                 << std::endl
		<< "   delK="  << deltaKE               << std::endl
		<< "   delE="  << delE                  << std::endl
		<< std::endl;
    }

    if (p2->dattrib[use_cons] - p2E > 0.0 and p2->dattrib[use_cons]>0.0) {
      std::cout << "P2 above zero: dif=" 
		<< p2->dattrib[use_cons] - p2E          << std::endl
		<< "    x_f="  << p2->dattrib[use_cons] << std::endl
		<< "    x_i="  << p2E                   << std::endl
		<< "    exs="  << Exs                   << std::endl
		<< "    misE=" << missE                 << std::endl
		<< "    delK=" << deltaKE               << std::endl
		<< "    delE=" << delE                  << std::endl
		<< std::endl;
    }

  }


  // Enforce electron equipartition
  //
  if (equiptn and use_elec) {

    // Electron from Particle 2
    //
    if (interFlag > 100 and interFlag < 200) {

      m1 = atomic_weights[Z2];
      m2 = atomic_weights[0 ];
      Mt = m1 + m2;
      Mu = m1 * m2 / Mt;

      KE1i = KE2i = 0.0;
      KE1f = KE2f = 0.0;

      double cost = 0.0, VC2 = 0.0, VR2 = 0.0;

      for (size_t k=0; k<3; k++) {
	KE1i += p2->vel[k] * p2->vel[k];
	KE2i += p2->dattrib[use_elec+k] * p2->dattrib[use_elec+k];
	cost += p2->vel[k] * p2->dattrib[use_elec+k];

	vcom[k] = (m1*p2->vel[k] + m2*p2->dattrib[use_elec+k])/Mt;
	vrel[k] = p2->vel[k] - p2->dattrib[use_elec+k];

	VC2    += vcom[k] * vcom[k];
	VR2    += vrel[k] * vrel[k];
      }

      if (KE1i > 0.0 and KE2i > 0.0) cost /= sqrt(KE1i * KE2i);

      double dmr   = cost / (m1 - m2);
      double gamma = 1.0 + 4.0*Mt*Mu*dmr*dmr;
      double E0    = 0.5*Mt*VC2 + 0.5*Mu*VR2;

      double gamP  = 1.0 + sqrt(1.0 - 1.0/gamma);
      double gamN  = 1.0 - sqrt(1.0 - 1.0/gamma);

      double virP  = 
	(VC2 - E0/Mt*gamN)*(VC2 - E0/Mt*gamN) +
	(VR2 - E0/Mu*gamP)*(VR2 - E0/Mu*gamP) ;
	
      double virN  = 
	(VC2 - E0/Mt*gamP)*(VC2 - E0/Mt*gamP) +
	(VR2 - E0/Mu*gamN)*(VR2 - E0/Mu*gamN) ;
	
      double vcfac = 0.0, vrfac = 0.0;

      if (virP > virN) {
	vcfac = sqrt(E0/Mt*gamN);
	vrfac = sqrt(E0/Mu*gamP);
      } else {
	vcfac = sqrt(E0/Mt*gamP);
	vrfac = sqrt(E0/Mu*gamN);
      }

      if (VC2>0.0) {
	for (size_t k=0; k<3; k++) vcom[k] /= sqrt(VC2);
      } else {
	double cos_th = 1.0 - 2.0*(*unit)();
	double sin_th = sqrt(1.0 - cos_th*cos_th);
	double phi    = 2.0*M_PI*(*unit)();
	vcom[0] = cos_th;
	vcom[1] = sin_th*cos(phi);
	vcom[2] = sin_th*sin(phi);
      }

      if (VR2>0.0) {
	for (size_t k=0; k<3; k++) vrel[k] /= sqrt(VR2);
      } else {
	double cos_th = 1.0 - 2.0*(*unit)();
	double sin_th = sqrt(1.0 - cos_th*cos_th);
	double phi    = 2.0*M_PI*(*unit)();
	vrel[0] = cos_th;
	vrel[1] = sin_th*cos(phi);
	vrel[2] = sin_th*sin(phi);
      }

      for (size_t k=0; k<3; k++) {
	p2->vel[k]              = vcom[k]*vcfac + m2/Mt * vrel[k]*vrfac;
	p2->dattrib[use_elec+k] = vcom[k]*vcfac - m1/Mt * vrel[k]*vrfac;

	KE1f += p2->vel[k] * p2->vel[k];
	KE2f += p2->dattrib[use_elec+k] * p2->dattrib[use_elec+k];
      }

      KE1i *= 0.5*m1;
      KE1f *= 0.5*m1;

      KE2i *= 0.5*m2;
      KE2f *= 0.5*m2;

      double KEi = KE1i + KE2i;
      double KEf = KE1f + KE2f;

      if ( fabs(KEi - KEf) > 1.0e-14*KEi ) {
	std::cout << "Test(1): keI=[" 
		  << std::setw(16) << KE1i << ", " 
		  << std::setw(16) << KE2i << "] keF=[" 
		  << std::setw(16) << KE1f << ", " 
		  << std::setw(16) << KE2f << "] vir=[" 
		  << std::setw(16) << virP << ", "
		  << std::setw(16) << virN << "] "
		  << std::endl;
      }

    } // end: electron from Particle 2

    // Electron from Particle 1
    //
    if (interFlag > 200 and interFlag < 300) {

      m1 = atomic_weights[Z1];
      m2 = atomic_weights[0 ];
      Mt = m1 + m2;
      Mu = m1 * m2 / Mt;

      KE1i = KE2i = 0.0;
      KE1f = KE2f = 0.0;

      double cost = 0.0, VC2 = 0.0, VR2 = 0.0;

      for (size_t k=0; k<3; k++) {
	KE1i += p1->vel[k] * p1->vel[k];
	KE2i += p1->dattrib[use_elec+k] * p1->dattrib[use_elec+k];
	cost += p1->vel[k] * p1->dattrib[use_elec+k];

	vcom[k] = (m1*p1->vel[k] + m2*p1->dattrib[use_elec+k])/Mt;
	vrel[k] = p1->vel[k] - p1->dattrib[use_elec+k];

	VC2    += vcom[k] * vcom[k];
	VR2    += vrel[k] * vrel[k];
      }

      if (KE1i > 0.0 and KE2i > 0.0) cost /= sqrt(KE1i * KE2i);

      double dmr   = cost / (m1 - m2);
      double gamma = 1.0 + 4.0*Mt*Mu*dmr*dmr;
      double E0    = 0.5*Mt*VC2 + 0.5*Mu*VR2;

      double gamP  = 1.0 + sqrt(1.0 - 1.0/gamma);
      double gamN  = 1.0 - sqrt(1.0 - 1.0/gamma);

      double virP  = 
	(VC2 - E0/Mt*gamN)*(VC2 - E0/Mt*gamN) +
	(VR2 - E0/Mu*gamP)*(VR2 - E0/Mu*gamP) ;
	
      double virN  = 
	(VC2 - E0/Mt*gamP)*(VC2 - E0/Mt*gamP) +
	(VR2 - E0/Mu*gamN)*(VR2 - E0/Mu*gamN) ;
	
      double vcfac = 0.0, vrfac = 0.0;

      if (virP > virN) {
	vcfac = sqrt(E0/Mt*gamN);
	vrfac = sqrt(E0/Mu*gamP);
      } else {
	vcfac = sqrt(E0/Mt*gamP);
	vrfac = sqrt(E0/Mu*gamN);
      }

      if (VC2>0.0) {
	for (size_t k=0; k<3; k++) vcom[k] /= sqrt(VC2);
      } else {
	double cos_th = 1.0 - 2.0*(*unit)();
	double sin_th = sqrt(1.0 - cos_th*cos_th);
	double phi    = 2.0*M_PI*(*unit)();
	vcom[0] = cos_th;
	vcom[1] = sin_th*cos(phi);
	vcom[2] = sin_th*sin(phi);
      }

      if (VR2>0.0) {
	for (size_t k=0; k<3; k++) vrel[k] /= sqrt(VR2);
      } else {
	double cos_th = 1.0 - 2.0*(*unit)();
	double sin_th = sqrt(1.0 - cos_th*cos_th);
	double phi    = 2.0*M_PI*(*unit)();
	vrel[0] = cos_th;
	vrel[1] = sin_th*cos(phi);
	vrel[2] = sin_th*sin(phi);
      }

      for (size_t k=0; k<3; k++) {
	p1->vel[k]              = vcom[k]*vcfac + m2/Mt * vrel[k]*vrfac;
	p1->dattrib[use_elec+k] = vcom[k]*vcfac - m1/Mt * vrel[k]*vrfac;

	KE1f += p1->vel[k] * p1->vel[k];
	KE2f += p1->dattrib[use_elec+k] * p1->dattrib[use_elec+k];
      }

      KE1i *= 0.5*m1;
      KE1f *= 0.5*m1;

      KE2i *= 0.5*m2;
      KE2f *= 0.5*m2;

      double KEi = KE1i + KE2i;
      double KEf = KE1f + KE2f;

      if ( fabs(KEi - KEf) > 1.0e-14*KEi ) {
	std::cout << "Test(1): keI=[" 
		  << std::setw(16) << KE1i << ", " 
		  << std::setw(16) << KE2i << "] keF=[" 
		  << std::setw(16) << KE1f << ", " 
		  << std::setw(16) << KE2f << "] vir=[" 
		  << std::setw(16) << virP << ", "
		  << std::setw(16) << virN << "] "
		  << std::endl;
      }

    }  // end: electron from Particle 2

  } // Equipartition stanza for electrons

  // Scatter electrons
  //
  if (esType == always and C1>1 and C2>1) {
    double vi = 0.0;
    for (int k=0; k<3; k++) {
      double d1 = p1->dattrib[use_elec+k];
      double d2 = p2->dattrib[use_elec+k];
      vcom[k] = 0.5*(d1 + d2);
      vi     += (d1 - d2) * (d1 - d2);
    }
    vi = sqrt(vi);

    double cos_th = 1.0 - 2.0*(*unit)();       // Cosine and sine of
    double sin_th = sqrt(1.0 - cos_th*cos_th); // Collision angle theta
    double phi    = 2.0*M_PI*(*unit)();	     // Collision angle phi
    
    vrel[0] = vi * cos_th;	  // Compute post-collision relative
    vrel[1] = vi * sin_th*cos(phi); // velocity for an elastic 
    vrel[2] = vi * sin_th*sin(phi); // interaction

    for (int k=0; k<3; k++) {
      p1->dattrib[use_elec+k] = vcom[k] + 0.5*vrel[k];
      p2->dattrib[use_elec+k] = vcom[k] - 0.5*vrel[k];
    }

    return 0;
  }

  // NOCOOL debugging
  //
  double fI1 = 0.0, fE1 = 0.0;
  double fI2 = 0.0, fE2 = 0.0;
  if (KE_DEBUG and NOCOOL) {
    for (size_t k=0; k<3; k++) {
      fI1 += p1->vel[k]*p1->vel[k];
      fI2 += p2->vel[k]*p2->vel[k];
      if (use_elec>=0) {
	if (C1>1) fE1 += p1->dattrib[use_elec+k]*p1->dattrib[use_elec+k];
	if (C2>1) fE2 += p2->dattrib[use_elec+k]*p2->dattrib[use_elec+k];
      }
    }
    fI1 *= 0.5*p1->mass;
    fI2 *= 0.5*p2->mass;
    fE1 *= 0.5*p1->mass*atomic_weights[0]/atomic_weights[Z1];
    fE2 *= 0.5*p2->mass*atomic_weights[0]/atomic_weights[Z2];

    double Einit = iI1 + iI2 + iE1 + iE2;
    double Efinl = fI1 + fI2 + fE1 + fE2;

    double testE = Einit - Efinl - NCXTRA;

    if (Z1==Z2)			// Add in energy loss/gain
      testE += Exs - delE - missE;
				// Correct for trace-algorithm excess
    else if ((C1==1 and C2==1) or electronic)
      testE -= deltaKE;

    if (fabs(testE) > DEBUG_THRESH*Einit )
      std::cout << "NC total ("<< m1 << "," << m2 << ") = " 
		<< std::setw(14) << testE
		<< ", flg=" << std::setw(6)  << interFlag
		<< ", dKE=" << std::setw(14) << Einit - Efinl
		<< ", com=" << std::setw(14) << deltaKE
		<< ", exs=" << std::setw(14) << Exs
		<< ", del=" << std::setw(14) << delE
		<< ", mis=" << std::setw(14) << missE
		<< ", NCX=" << std::setw(14) << NCXTRA
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

  //
  // Cross section scale factor
  //
  double scaleCrossSection = 0.0;

  for (auto sp : sCross[id]) {

    for (auto v : sCross[id][sp.first]) {

      // Mass fractions
      //
      double w1 = p1->dattrib[SpList[sp.first]];
      double w2 = p2->dattrib[SpList[sp.first]];

      scaleCrossSection += v * (w1 + w2) / atomic_weights[sp.first.first];
    }
  }

  // Test
  //
  scaleCrossSection /= csections[id][defaultKey][defaultKey] *
    1e-14 / (UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
  
  //-------------------------
  // VERBOSE DEBUG TEST
  //-------------------------
  //
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
      double N1 = Pr1 * Fn * p1->mass * scaleCrossSection;
      double N2 = Pr2 * Fn * p2->mass * scaleCrossSection;

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

	if (RECOMB_IP) {
	  double Xi = ch.IonList[lQ(Z, C-1)]->ip;
	  delE1 = Xi * N1;
	  delE2 = Xi * N2;
	}
	
	delE1 += kEe2[id] * N1;
	delE2 += kEe1[id] * N2;

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

  // Artifically prevent cooling by setting the energy removed from
  // the COM frame to zero
  //
  if (NOCOOL) delE = 0.0;

  // Convert energy loss to system units
  //
  delE = delE/UserTreeDSMC::Eunit;
  
  // Assign interaction energy variables
  //
  double totE, kEe = kEe1[id];

  // Diagnostic accumulation
  //
  totE = kE;			// KE
  
  // Warn if energy lost is greater than total energy available to
  // lose
  //
  if (frost_warning && delE > totE)
    std::cout << "delE > KE!! (" << delE << " > " << totE
	      << "), kEe  = "  << kEe
	      << " delE = " << delE/(eV*Mu*UserTreeDSMC::Munit*amu)
	      << std::endl;
  
  // Cooling rate diagnostic histogram
  //
  if (TSDIAG && delE>0.0) {
				// Histogram index
    int indx = (int)floor(log(totE/delE)/(log(2.0)*TSPOW) + 5);
				// Floor and ceiling
    if (indx<0 ) indx = 0;
    if (indx>10) indx = 10;
				// Add entry
    EoverT[id][indx] += Mt;
  }
  
  // Time step "cooling" diagnostic
  //
  if (use_delt>=0 && delE>0.0 && totE>0.0) {
    double dtE = totE/delE * spTau[id];
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
    if (totE > delE) {
      
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

      (*cr)          = 0.0;
      ret            = 1;	// Set error flag
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

void CollideIon::eEdbg()
{
  if (fabs(tot2[0] - tot2[1])/tot2[0] > 1.0e-14) {

    std::ofstream out(tpaths[itp++ % 4].c_str());
    for (auto m : data[0]) {
      out << std::setw(10) << m.first
	  << std::setw(14) << std::get<0>(m.second)
	  << std::setw(14) << std::get<1>(m.second)
	  << std::setw(14) << std::get<0>(data[1][m.first])
	  << std::setw(14) << std::get<1>(data[1][m.first])
	  << std::endl;
    }
    out << std::setw(10) << "***"
	<< std::setw(14) << tot1[0]
	<< std::setw(14) << tot2[0]
	<< std::setw(14) << tot1[1]
	<< std::setw(14) << tot2[1]
	<< std::endl;
  }
}

double CollideIon::electronEnergy(pCell* const cell, int dbg)
{
  double Eengy = 0.0;
  for (auto b : cell->bods) {
    Particle *p = c0->Tree()->Body(b);
    KeyConvert k(p->iattrib[use_key]);
    if (k.C() - 1 > 0) {
      double numb = p->mass/atomic_weights[k.Z()];
      for (unsigned j=0; j<3; j++) {
	double v = p->dattrib[use_elec+j];
	Eengy += 0.5 * v*v * numb;
      }
    }
  }

  if (dbg>=0) {

    data[dbg].clear();
    tot1[dbg] = 0.0;
    tot2[dbg] = 0.0;

    for (auto b : cell->bods) {
      Particle *p = c0->Tree()->Body(b);
      KeyConvert k(p->iattrib[use_key]);
      if (k.C() - 1 > 0) {
	double numb = p->mass/atomic_weights[k.Z()];
	double E = 0.0;
	for (unsigned j=0; j<3; j++) {
	  double v = p->dattrib[use_elec+j];
	  E += 0.5 * v*v * numb;
	}
	std::get<0>(data[dbg][b]) = 2.0*E/numb;
	std::get<1>(data[dbg][b]) = E;
	tot1[dbg] += 2.0*E/numb;
	tot2[dbg] += E;
      }
    }
  }

  return Eengy * atomic_weights[0];
}

void CollideIon::finalize_cell(pHOT* const tree, pCell* const cell, 
			       sKeyDmap* const Fn, double kedsp, double tau,
			       int id)
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
  } // end: Trace
  
  
  // Do electron interactions separately
  //
  if ( (aType == Direct or aType == Weight) and use_elec>=0 and
       esType != always and esType != none) {

    const double cunit = 1e-14/(UserTreeDSMC::Lunit*UserTreeDSMC::Lunit);
    std::vector<unsigned long> bods;
    double eta = 0.0, crsvel = 0.0;
    double volc = cell->Volume();
    double me   = atomic_weights[0]*amu;

    // Momentum diagnostic distribution
    //
    std::vector<double> pdif;

    // Compute list of particles in cell with electrons
    //
    for (auto i : cell->bods) {
      Particle *p = cell->Body(i);
      KeyConvert k(p->iattrib[use_key]);
      if (k.C()>1) {
	bods.push_back(i);
	eta += ZMList[k.Z()]/atomic_weights[k.Z()] * 
	  (*Fn)[k.getKey()] * (k.C()-1);
      }
    }
    
    // Sample cell
    //
    pCell *samp = cell->sample;

    if (samp)
      crsvel = ntcdb[samp->mykey].CrsVel(electronKey, 0.95);
    else
      crsvel = ntcdb[cell->mykey].CrsVel(electronKey, 0.95);
    
    // Probability of an interaction of between particles of type 1
    // and 2 for a given particle of type 2
    //
    double    Prob = eta * cunit * crsvel * tau / volc;
    size_t   nbods = bods.size();
    double   selcM = 0.5 * nbods * (nbods-1) *  Prob;
    unsigned nselM = static_cast<unsigned>(floor(selcM+0.5));

    if (esType == limited or esType == fixed) 
      nselM = std::min<unsigned>(nselM, esNum);

    for (unsigned n=0; n<nselM; n++) {
    
      // Pick two particles with electrons at random out of this cell.
      //
      size_t l1, l2;

      l1 = static_cast<size_t>(floor((*unit)()*nbods));
      l1 = std::min<size_t>(l1, nbods-1);

      l2 = static_cast<size_t>(floor((*unit)()*(nbods-1)));
      l2 = std::min<size_t>(l2, nbods-2);
	  
      if (l2 >= l1) l2++;

      // Get index from body map for the cell
      //
      Particle* p1 = cell->Body(bods[l1]);
      Particle* p2 = cell->Body(bods[l2]);

      KeyConvert k1(p1->iattrib[use_key]);
      KeyConvert k2(p2->iattrib[use_key]);

      if (SAME_ELEC_SCAT) if (k1.Z() != k2.Z()) continue;

      // Swap particles so that p2 is the trace element
      //
      if (p1->mass/atomic_weights[k1.Z()] < 
	  p2->mass/atomic_weights[k2.Z()]) 
	{
	  // Swap the particle pointers
	  //
	  Particle *pT = p1;
	  p1 = p2;
	  p2 = pT;

	  // Reassign the keys and species indices
	  //
	  k1 = KeyConvert(p1->iattrib[use_key]);
	  k2 = KeyConvert(p2->iattrib[use_key]);
	}

      double ne1 = k1.C() - 1;
      double ne2 = k2.C() - 1;

      // Find the trace ratio
      //
      double Wa = p1->mass / atomic_weights[k1.Z()];
      double Wb = p2->mass / atomic_weights[k2.Z()];
      double  q = Wb / Wa;

      double ma = atomic_weights[0];
      double mb = atomic_weights[0];
      double mt = ma + mb;

      // Calculate pair's relative speed (pre-collision)
      //
      vector<double> vcom(3), vrel(3), v1(3), v2(3);
      for (int k=0; k<3; k++) {
	v1[k] = p1->dattrib[use_elec+k];
	v2[k] = p2->dattrib[use_elec+k];
	vcom[k] = (ma*v1[k] + mb*v2[k])/mt;
	vrel[k] = v1[k] - v2[k];
      }
      
      // Compute relative speed
      //
      double vi = 0.0;
      for (auto v : vrel) vi += v*v;

      // No point in inelastic collsion for zero velocity . . . 
      //
      if (vi == 0.0) continue;

      // Relative velocity 
      //
      vi = sqrt(vi);

      double cr = vi * UserTreeDSMC::Vunit;

      // Kinetic energy in eV
      //
      double kEee = 0.25 * me * cr * cr / eV;

      // Compute the cross section
      //
      double scrs = 0.0;
      
      // Mean interparticle spacing
      // 
      double ips = pow(volc/nbods, 0.333333) * UserTreeDSMC::Lunit * 1.0e7;
      
      // Collision flag
      //
      bool ok = true;
      
      if (esType == classical or esType == limited) {

	if (use_elec >=0 and ne1 > 0 and ne2 > 0) {
	  double b = 0.5*esu*esu /
	    std::max<double>(kEee*eV, FloorEv*eV) * 1.0e7; // nm
	  b = std::min<double>(b, ips);
	  scrs = M_PI*b*b * ne1 * ne2 * logL;
	}

	// Accept or reject candidate pair according to relative speed
	//
	double prod = vi * scrs;
	double targ = ntcdb[samp->mykey].Prob(electronKey, prod);

	ok = (targ > (*unit)() );
	
	// Over NTC max average
	//
	if (targ >= 1.0) elecOvr[id]++;

	// Used / Total
	//
	if (ok)         elecAcc[id]++; elecTot[id]++;

	// Update v_max and cross_max for NTC
	//
#pragma omp critical
	ntcdb[samp->mykey].Add(electronKey, prod);
      }


      double deltaKE = 0.0, dKE = 0.0;
      double KEi1 = 0.0, KEi2 = 0.0; // For debugging

      if (KE_DEBUG) {

	for (auto v : v1) KEi1 += v*v;
	for (auto v : v2) KEi2 += v*v;

	KEi1 *= 0.5*Wa*ma;
	KEi2 *= 0.5*Wb*mb;
      }

      // Scatter
      //
      if (ok) {
	double cos_th = 1.0 - 2.0*(*unit)();       // Cosine and sine of
	double sin_th = sqrt(1.0 - cos_th*cos_th); // Collision angle theta
	double phi    = 2.0*M_PI*(*unit)();        // Collision angle phi
    
	vrel[0] = vi * cos_th;	        // Compute post-collision
	vrel[1] = vi * sin_th*cos(phi);	// relative velocity for an
	vrel[2] = vi * sin_th*sin(phi);	// elastic interaction

	// Explicit energy conservation using splitting
	//
	if (ExactE) {

	  double vrat = 1.0;

	  std::vector<double> uu(3), vv(3), w1(v1);
	  for (size_t k=0; k<3; k++) {
				// New velocities in COM
	    uu[k] = vcom[k] + 0.5*vrel[k];
	    vv[k] = vcom[k] - 0.5*vrel[k];
	  }

	  if (q < 1.0) {

	    double v1i2 = 0.0, b1f2 = 0.0, v2i2 = 0.0, b2f2 = 0.0;
	    for (size_t k=0; k<3; k++) {
				// Compute energies and angle
	      v1i2 += v1[k]*v1[k];
	      v2i2 += v2[k]*v2[k];
	      b1f2 += uu[k]*uu[k];
	      b2f2 += vv[k]*vv[k];
	    }

	    if (AlgOrth) {
	      // Cross product to determine orthgonal direction
	      //
	      w1[0] = uu[1]*v1[2] - uu[2]*v1[1];
	      w1[1] = uu[2]*v1[0] - uu[0]*v1[2];
	      w1[2] = uu[0]*v1[1] - uu[1]*v1[0];

	      // Normalize
	      //
	      double wnrm = 0.0;
	      for (auto   v : w1) wnrm += v*v;
	      for (auto & v : w1) v *= 1.0/sqrt(wnrm);
	      
	      vrat  = sqrt((1.0 - q)*(q*b1f2 + v1i2))/(1.0 - q);
	      
	    } else {

	      double qT = 0.0;
	      for (size_t k=0; k<3; k++) qT += v1[k]*uu[k];
	    
	      if (v1i2 > 0.0 and b1f2 > 0.0) qT *= q/sqrt(v1i2 * b1f2);
	  
	      double vh1f  = 
		( -sqrt(b1f2)*qT + sqrt(qT*qT*b1f2 + (1.0 - q)*(q*b1f2 + v1i2)) )/(1.0 - q);
	    
	      vrat = vh1f / sqrt(v1i2);
	    }

	  }

	  // New velocities in inertial frame
	  //
	  std::vector<double> u1(3), u2(3);
	  for (size_t k=0; k<3; k++) {
	    u1[k] = (1.0 - q)*w1[k]*vrat + q*uu[k];
	    u2[k] = vv[k];
	  }

	  // These are all for diagnostics
	  //
	  double pi2 = 0.0, dp2 = 0.0, Ebeg = 0.0, Efin = 0.0;

	  // DIAGNOSTIC: energies
	  for (int k=0; k<3; k++) {
	    Ebeg += 0.5*Wa*ma*v1[k]*v1[k] + 0.5*Wb*mb*v2[k]*v2[k];
	    Efin += 0.5*Wa*ma*u1[k]*u1[k] + 0.5*Wb*mb*u2[k]*u2[k];
	  }

	  // Assign new electron velocities
	  //
	  for (int k=0; k<3; k++) {
	    p1->dattrib[use_elec+k] = u1[k];
	    p2->dattrib[use_elec+k] = u2[k];

	    // DIAGNOSTIC: initial and final momenta
	    double dpi = Wa*ma*v1[k] + Wb*mb*v2[k];
	    double dpf = Wa*ma*u1[k] + Wb*mb*u2[k];

	    // DIAGNOSTIC: rms momentum difference
	    dp2  += (dpi - dpf)*(dpi - dpf); 
	    pi2  += dpi*dpi;		     
	  }

	  // Check for energy conservation
	  //
	  if (DebugE) momD[id].push_back(sqrt(dp2/pi2));

	  if ( fabs(Efin - Ebeg) > 1.0e-12*(Ebeg) ) {
	    std::cout << "Broken energy conservation,"
		      << " Ebeg="  << Ebeg
		      << " Efin="  << Efin
		      << " Edif="  << Efin/Ebeg - 1.0
		      << " pcons=" << sqrt(dp2/pi2)
		      << std::endl;
	  }

	  // Upscale electron energy
	  //
	  if (TRACE_ELEC and k1.Z() == k2.Z()) {
	    double delE = p1->dattrib[use_elec+3] + p2->dattrib[use_elec+3];
	    if (delE > 0.0) {
	      p1->dattrib[use_elec+3] = p2->dattrib[use_elec+3] = 0.0;

	      double m1    = p1->mass*atomic_weights[0]/atomic_weights[k1.Z()];
	      double m2    = p2->mass*atomic_weights[0]/atomic_weights[k2.Z()];
	      double mt    = m1 + m2;
	      double mu    = m1 * m2 / mt;
	      double KEcom = 0.0;

	      for (int k=0; k<3; k++) {
		vcom[k] = (m1*u1[k] + m2*u2[k]);
		vrel[k] = u1[k] - u2[k];
		KEcom  += vrel[k] * vrel[k];
	      }
	      KEcom *= 0.5*mu;

	      if (KEcom>0.0) {
		double vfac = sqrt(1.0 + delE/KEcom);
		for (int k=0; k<3; k++) {
		  p1->dattrib[use_elec+k] = vcom[k] + m2/mt*vrel[k]*vfac;
		  p2->dattrib[use_elec+k] = vcom[k] - m1/mt*vrel[k]*vfac;
		}
	      }
	      
	    } // end: delE > 0.0
	    
	  } // end: TRACE_ELEC

	} // end: ExactE
	
	// Explicit momentum conservation
	//
	else {

	  bool equal = fabs(q - 1.0) < 1.0e-14;

	  double vfac = 1.0;
	  if (equal) {
	    double KE0 = 0.5*Wa*ma*mb/mt*vi*vi;
	    dKE = p1->dattrib[use_elec+3] + p2->dattrib[use_elec+3];
	    vfac = sqrt(1.0 + dKE/KE0);
	    p1->dattrib[use_elec+3] = p2->dattrib[use_elec+3] = 0.0;
	  }

	  double qKEfac = 0.5*Wa*ma*q*(1.0 - q);
	  for (int k=0; k<3; k++) {
	    double v0 = vcom[k] + mb/mt*vrel[k]*vfac;
	    deltaKE += (v0 - v1[k])*(v0 - v1[k]) * qKEfac;
	    p1->dattrib[use_elec+k] = (1.0 - q)*v1[k] + q*v0;
	    p2->dattrib[use_elec+k] = vcom[k] - ma/mt*vrel[k]*vfac;
	  }
				// Correct energy for conservation
	  if (!equal) p1->dattrib[use_elec+3] += deltaKE;
	}
      }

      if (KE_DEBUG) {

	double KEf1 = 0.0;
	double KEf2 = 0.0;

	for (int k=0; k<3; k++) {
	  double v1 = p1->dattrib[use_elec+k];
	  double v2 = p2->dattrib[use_elec+k];
	  KEf1 += v1*v1;
	  KEf2 += v2*v2;
	}

	KEf1 *= 0.5*Wa*ma;
	KEf2 *= 0.5*Wb*mb;
	
	double KEi = KEi1 + KEi2;
	double KEf = KEf1 + KEf2;

	double testE = KEi - KEf - deltaKE + dKE;
	
	if (fabs(testE) > DEBUG_THRESH*KEi) {
	  std::cout << std::endl << std::string(70, '-') << std::endl
		    << "Total elec ("
		    << k1.Z() << "," 
		    << k2.Z() << ") = "
		    << std::setw(14) << testE
		    << ", KEi=" << std::setw(14) << KEi
		    << ", KEf=" << std::setw(14) << KEf
		    << ", exs=" << std::setw(14) << deltaKE
		    << ", dKE=" << std::setw(14) << dKE
		    << std::endl;
	}

      }

    } // loop over particles


  } // end: Direct or Weight for use_elec>=0


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
	    << "# N(nn)         number of neut-neut scat " << std::endl
	    << "# N(ne)         number of neut-elec scat " << std::endl
	    << "# N(ie)         number of ion-elec scat  " << std::endl
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
	  size_t w =16*12, l = sout.str().size();
	  sout2 << std::setw((w-l)/2) << ' ' << sout.str();
	  out   << std::setw(w) << sout2.str() << " | ";
	}
	out << std::setw(2*12) << ' ' << " |" << std::endl;

				// Header line
	out << std::setfill('-') << std::right;
	out << "#" << std::setw(11+12) << '+' << " | ";
	for (auto it : *this) {
	  for (int i=0; i<16; i++) out << std::setw(12) << '+';
	  out << " | ";
	}
	out << std::setw(12) << '+' << std::setw(12) << '-' << " |"
	    << std::setfill(' ') << std::endl;

				// Column labels
	out << "#" 
	    << std::setw(11) << "Time |"
	    << std::setw(12) << "Temp |" << " | ";
	for (auto it : *this) {
	  out << std::setw(12) << "N(nn) |"
	      << std::setw(12) << "N(ne) |"
	      << std::setw(12) << "N(ie) |"
	      << std::setw(12) << "N(ff) |"
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
	  for (size_t l=0; l<16; l++) {
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
	  for (int i=0; i<16; i++) out << std::setw(12) << '+';
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
	out << std::setw(12) << std::get<0>(ctd->nn_s)
	    << std::setw(12) << std::get<0>(ctd->ne_s)
	    << std::setw(12) << std::get<0>(ctd->ie_s)
	    << std::setw(12) << std::get<0>(ctd->ff_s)
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

	in.getline(line, nline);

	if (in.good()) {
	  std::istringstream sz(line);
	  sz >> use_elec;
	} else {
	  nOK = 1;		// Can't read position flag
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
  sKeyDmap            eta, densM, densN, collP, nsigmaM, ncrossM;
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
  
  for (auto it1 : c->count) {

    // Only compute if particles of this species is in the cell
    //
    if (it1.second) {
      speciesKey i1 = it1.first;
      
      // Trace weight, Eta_b
      //
      eta[i1] = ZMList[i1.first] / atomic_weights[i1.first];

      // Mass density scaled by atomic weight in amu.  In the
      // algorithm notes, this is N_b * Eta_b / V.
      //
      densM[i1] = c->Mass(i1) / atomic_weights[i1.first] / volc;
      
      // Number density of superparticles
      //
      densN[i1] = static_cast<double>(c->Count(i1))/volc;
    }
  }

  if (DEBUG_SL) {

    std::cout << std::endl
	      << std::string(70, '-')     << std::endl
	      << "Cell stats"
	      << ", #=" << c->bods.size() << std::endl
	      << std::string(70, '-')     << std::endl
	      << std::setw(10) << "Species"
	      << std::setw(16) << "eta"
	      << std::setw(16) << "n dens"
	      << std::setw(16) << "m dens"
	      << std::setw(16) << "sp mass"
	      << std::setw(10) << "n count"
	      << std::endl
	      << std::setw(10) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(16) << "---------"
	      << std::setw(10) << "---------"
	      << std::endl;

    for (auto it : c->count) {
      std::ostringstream sout;
      sout << "(" << it.first.first << ", " << it.first.second << ")";
      std::cout << std::setw(10) << sout.str()
		<< std::setw(16) << eta  [it.first]
		<< std::setw(16) << densN[it.first]
		<< std::setw(16) << densM[it.first]
		<< std::setw(16) << c->Mass(it.first)
		<< std::setw(10) << c->Count(it.first)
		<< std::endl;
    }

    std::cout << std::endl
	      << std::string(70, '-')    << std::endl
	      << "Interaction stats"
	      << ", eVel=" << Evel[id] 
	      << ", crm="  << crm        << std::endl
	      << std::string(70, '-')    << std::endl
	      << std::setw(20) << "Species"
	      << std::setw(16) << "Cross"
	      << std::endl
	      << std::setw(20) << "---------"
	      << std::setw(16) << "---------"
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
	std::cout << std::setw(20) << sout.str()
		  << std::setw(16) << csections[id][k1][k2] / cunit
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
	    crossT      *= (*Fn)[i2] * eta[i2];
	    ncrossM[i1] += crossT;
	    nsigmaM[i1] += densN[i2] * crossT;
	  } else {
	    crossT      *= (*Fn)[i1] * eta[i1];
	    ncrossM[i2] += crossT;
	    nsigmaM[i2] += densN[i1] * crossT;
	  }
      
	  // So, ncrossM is the superparticle cross section for each species

	  // Sanity check debugging
	  //
	  if (csections[id][i1][i2] <= 0.0 || std::isnan(csections[id][i1][i2])) {
	    cout << "INVALID CROSS SECTION! :: " << csections[id][i1][i2]
		 << " #1 = (" << i1.first << ", " << i1.second << ")"
		 << " #2 = (" << i2.first << ", " << i2.second << ")";

	    csections[id][i1][i2] = 0.0; // Zero out
	  }
	  
	  // Sanity check debugging
	  //
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

      // Sanity check debugging
      //
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
	    crsvel = ntcdb[samp->mykey].CrsVel(k, 0.95);
	  else
	    crsvel = ntcdb[c->mykey].CrsVel(k, 0.95);
	  
	  // Probability of an interaction of between particles of type 1
	  // and 2 for a given particle of type 2
	  //
	  double Prob = 0.0;

	  if (densM[i1]>=densM[i2]) {
	    Prob = (*Fn)[i2] * eta[i2] * cunit * crsvel * tau / volc;
	  } else {
	    Prob = (*Fn)[i1] * eta[i1] * cunit * crsvel * tau / volc;
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
		cv1 = ntcdb[samp->mykey].CrsVel(k, 0.50);
		cv2 = ntcdb[samp->mykey].CrsVel(k, 0.90);
		cv3 = ntcdb[samp->mykey].CrsVel(k, 0.95);
	      } else {
		cv1 = ntcdb[c->mykey].CrsVel(k, 0.50);
		cv2 = ntcdb[c->mykey].CrsVel(k, 0.90);
		cv3 = ntcdb[c->mykey].CrsVel(k, 0.95);
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
	      crsvel = ntcdb[samp->mykey].CrsVel(k, 0.95);
	    else
	      crsvel = ntcdb[c->mykey].CrsVel(k, 0.95);

	    double Prob0 = 0.0, Prob1 = 0.0;

	    if (densM[i1]>=densM[i2]) {
	      Prob0 = densM[i2] * (*Fn)[i2] * cunit * crsvel * tau;
	      Prob1 = nsigmaM[i2] * crm * tau;
	    } else {
	      Prob0 = densM[i1] * (*Fn)[i1] * cunit * crsvel * tau;
	      Prob1 = nsigmaM[i1] * crm * tau;
	    }
	    
	    std::cout << "(" 
		      << std::setw(2)  << i1.first << ","
		      << std::setw(2)  << i1.second << ") ("
		      << std::setw(2)  << i2.first << ","
		      << std::setw(2)  << i2.second << ")  "
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
	      << "  MFP/L = "       << meanLambda/pow(volc, 0.333333333)
	      << "  totalNsel = "   << totalNsel
	      << std::endl << std::endl;
  }
  
  return nselM;
}

sKey2Umap CollideIon::generateSelectionTrace
(pCell* const c, sKeyDmap* const Fn, double crm, double tau, int id,
 double& meanLambda, double& meanCollP, double& totalNsel)
{
  speciesKey key(defaultKey);
    
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

    double mass  = 0.0;

    consE = 0.0;
    consG = 0.0;
    totlE = 0.0;
    tempM = 0.0;
    tempE = 0.0;
    elecE = 0.0;
    
    specI.clear();
    specE.clear();

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

      if (std::isnan(KEtot)) {	// For debugging, obviously
	std::cout << "NaN in sample cell stats" << std::endl;
	cell->sample->KE(KEtot, KEdsp);
      }

      if (use_cons >= 0) {
	for (auto b : cell->bods) {
	  consE += c0->Tree()->Body(b)->dattrib[use_cons];
	  if (use_elec>=0)
	    consG += c0->Tree()->Body(b)->dattrib[use_elec+3];
	}
      }
      
      cType::iterator ft = ETcache.find(cell->sample->mykey);

      if (use_elec >= 0) {
	
	if (ft != ETcache.end()) {
	  T = ft->second;
	} else {
	  typedef std::tuple<double, double> dtup;
	  const dtup zero(0.0, 0.0);
	  std::vector<dtup> vel(3, zero);

	  double count = 0.0;
	  for (auto c : cell->sample->children) {
	    for (auto b : c.second->bods) {
	      Particle *p = c0->Tree()->Body(b);
	      KeyConvert k(p->iattrib[use_key]);
	      unsigned short ne = k.C() - 1;
	      double numb = p->mass/atomic_weights[k.Z()];
	      if (ne) {
		for (unsigned k=0; k<3; k++) {
		  double v = p->dattrib[use_elec+k];
		  std::get<0>(vel[k]) += v   * numb;
		  std::get<1>(vel[k]) += v*v * numb;
		}
		count += numb;
	      }
	    }
	  }
	    
	  double dispr = 0.0;
	    
	  // Temp computation
	  // ----------------
	  // number of atoms      = w_i  = m_i/(mu_i * m_a)
	  // elec mean velocity   = sv1  = sum_i (w_i v)
	  // elec mean vel^2      = sv2  = sum_i (w_i v^2)
	  // summed number        = sn   = sum_i (w_i)
	  // elec specific KE     = disp = sv2 - sv1*sv1/sn
	  // total elecron KE     = KE   = mu_e*m_a * disp
	  // number of elecrons   = N    = sum_i (m_i/(mu_i * m_a)) = sn
	  // total elecron KE     =        3/2 N k T
	  // KE prefactor         = Tfac = 2 * m_a/(3*k)
	  // ----------------
	  // Solve for T
	  // ----------------
	  // T = [2/(3*k) * m_a] * KE/sn
	  //   = [2/(3*k) * m_a] * mu_e * disp / sn
	  //   = Tfac * mu_e * disp / sn
	  
	  if (count > 0.0) {
	    for (auto v : vel) {
	      double v1 = std::get<0>(v);
	      dispr += 0.5*(std::get<1>(v) - v1*v1/count);
	    }
	    T = ETcache[cell->sample->mykey] = 
	      Tfac * atomic_weights[0] * dispr/count;
	  } else {
	    T = ETcache[cell->sample->mykey] = 0.0;
	  }
	}

				// Mass-weighted temperature
	tempE += cell->Mass() * T; 

	// Compute total electron energy in this cell
	//
	elecE += electronEnergy(cell);

	// Compute electron energy per element
	//
	for (auto b : cell->bods) {
	  Particle *p = c0->Tree()->Body(b);
	  unsigned Z  = KeyConvert(p->iattrib[use_key]).Z();
	  double num  = p->mass / atomic_weights[Z];
	  
	  if (specI.find(Z) == specI.end()) specI[Z] = ZTup(0, 0);
	  double v2 = 0.0;
	  for (auto v : p->vel) v2 += v*v;
	  std::get<0>(specI[Z]) += 0.5*num*atomic_weights[Z]*v2;
	  std::get<1>(specI[Z]) += num;
	  
	  if (KeyConvert(p->iattrib[use_key]).C()==1) continue;
	  if (specE.find(Z) == specE.end()) specE[Z] = ZTup(0, 0);
	  v2  = 0.0;
	  for (int j=0; j<3; j++) {
	    double v = p->dattrib[use_elec+j];
	    v2 += v*v;
	  }
	  std::get<0>(specE[Z]) += 0.5*num*atomic_weights[0]*v2;
	  std::get<1>(specE[Z]) += num;
	}

      } // end: use_elec>=0

    } // end: cell loop
    

    // Send values to root
    //
    double val1, val2, val3 = 0.0, val4 = 0.0, val5 = 0.0;
    double val6 = 0.0, val7 = 0.0;
    
    if (COLL_SPECIES) {
      for (int t=1; t<nthrds; t++) {
	for (auto s : collCount[t]) {
	  if (collCount[0].find(s.first) == collCount[0].end()) 
	    collCount[0][s.first]  = s.second;
	  else {
	    collCount[0][s.first][0] += s.second[0];
	    collCount[0][s.first][1] += s.second[1];
	  }
	}
      }
    }

    for (int i=1; i<numprocs; i++) {

      if (i == myid) {
				// Mass
	MPI_Send(&mass,  1, MPI_DOUBLE, 0, 331, MPI_COMM_WORLD);
				// Temp
	MPI_Send(&tempM, 1, MPI_DOUBLE, 0, 332, MPI_COMM_WORLD);

	MPI_Send(&consE, 1, MPI_DOUBLE, 0, 333, MPI_COMM_WORLD);
	MPI_Send(&consG, 1, MPI_DOUBLE, 0, 334, MPI_COMM_WORLD);
	MPI_Send(&totlE, 1, MPI_DOUBLE, 0, 335, MPI_COMM_WORLD);

				// Energies
	if (use_elec >= 0) {
	  MPI_Send(&tempE, 1, MPI_DOUBLE, 0, 336, MPI_COMM_WORLD);
	  MPI_Send(&elecE, 1, MPI_DOUBLE, 0, 337, MPI_COMM_WORLD);

				// Local ion map size
	  int sizm = specE.size();
	  MPI_Send(&sizm,  1, MPI_INT,    0, 338, MPI_COMM_WORLD);

				// Send local ion map
	  for (auto i : specI) {
	    unsigned short Z = i.first;
	    MPI_Send(&Z, 1, MPI_UNSIGNED_SHORT, 0, 339, MPI_COMM_WORLD);
	    double E = std::get<0>(i.second);
	    double N = std::get<1>(i.second);
	    MPI_Send(&E, 1, MPI_DOUBLE,         0, 340, MPI_COMM_WORLD);
	    MPI_Send(&N, 1, MPI_DOUBLE,         0, 341, MPI_COMM_WORLD);
	  }

				// Local electron map size
	  sizm = specE.size();
	  MPI_Send(&sizm,  1, MPI_INT,    0, 342, MPI_COMM_WORLD);

				// Send local electron map
	  for (auto e : specE) {
	    unsigned short Z = e.first;
	    MPI_Send(&Z, 1, MPI_UNSIGNED_SHORT, 0, 343, MPI_COMM_WORLD);
	    double E = std::get<0>(e.second);
	    double N = std::get<1>(e.second);
	    MPI_Send(&E, 1, MPI_DOUBLE,         0, 344, MPI_COMM_WORLD);
	    MPI_Send(&N, 1, MPI_DOUBLE,         0, 345, MPI_COMM_WORLD);
	  }

	} // end: use_elec>=0

	if (COLL_SPECIES) {
	  for (auto s : collCount[0]) {
	    speciesKey k1 = s.first.first;
	    speciesKey k2 = s.first.second;
	    MPI_Send(&k1.first,    1, MPI_UNSIGNED_SHORT, 0, 346, 
		     MPI_COMM_WORLD);
	    MPI_Send(&k1.second,   1, MPI_UNSIGNED_SHORT, 0, 347, 
		     MPI_COMM_WORLD);
	    MPI_Send(&k2.first,    1, MPI_UNSIGNED_SHORT, 0, 348, 
		     MPI_COMM_WORLD);
	    MPI_Send(&k2.second,   1, MPI_UNSIGNED_SHORT, 0, 349, 
		     MPI_COMM_WORLD);
	    MPI_Send(&s.second[0], 1, MPI_UNSIGNED_LONG,  0, 350,
		       MPI_COMM_WORLD);
	    MPI_Send(&s.second[1], 1, MPI_UNSIGNED_LONG,  0, 351,
		       MPI_COMM_WORLD);
	  }
	  unsigned short Z = 255;
	  MPI_Send(&Z, 1, MPI_UNSIGNED_SHORT, 0, 346, MPI_COMM_WORLD);
	}

      }	// end: myid>0

				// Root receives from Node i
      if (0 == myid) {

	MPI_Recv(&val1, 1, MPI_DOUBLE, i, 331, MPI_COMM_WORLD, 
		 MPI_STATUS_IGNORE);
	MPI_Recv(&val2, 1, MPI_DOUBLE, i, 332, MPI_COMM_WORLD, 
		 MPI_STATUS_IGNORE);

	MPI_Recv(&val3, 1, MPI_DOUBLE, i, 333, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	MPI_Recv(&val4, 1, MPI_DOUBLE, i, 334, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	MPI_Recv(&val5, 1, MPI_DOUBLE, i, 335, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);

	if (use_elec >= 0) {
	  MPI_Recv(&val6, 1, MPI_DOUBLE, i, 336, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  MPI_Recv(&val7, 1, MPI_DOUBLE, i, 337, MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
      
	  int sizm;
				// Receive ion map size
	  MPI_Recv(&sizm, 1, MPI_INT, i, 338, MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  
				// Receive ion map
	  for (int j=0; j<sizm; j++) {

	    unsigned short Z;
	    double   E;
	    double   N;

	    MPI_Recv(&Z, 1, MPI_UNSIGNED_SHORT, i, 339, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);
	    MPI_Recv(&E, 1, MPI_DOUBLE,         i, 340, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);
	    MPI_Recv(&N, 1, MPI_DOUBLE,         i, 341, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);

	    if (specI.find(Z) == specE.end()) specI[Z] = ZTup(0, 0);
	    std::get<0>(specI[Z]) += E;
	    std::get<1>(specI[Z]) += N;
	  }

				// Receive electron map size
	  MPI_Recv(&sizm, 1, MPI_INT, i, 342, MPI_COMM_WORLD, 
		   MPI_STATUS_IGNORE);
	  
				// Receive electron map
	  for (int j=0; j<sizm; j++) {

	    unsigned short Z;
	    double   E;
	    double   N;

	    MPI_Recv(&Z, 1, MPI_UNSIGNED_SHORT, i, 343, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);
	    MPI_Recv(&E, 1, MPI_DOUBLE,         i, 344, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);
	    MPI_Recv(&N, 1, MPI_DOUBLE,         i, 345, MPI_COMM_WORLD,
		     MPI_STATUS_IGNORE);

	    if (specE.find(Z) == specE.end()) specE[Z] = ZTup(0, 0);
	    std::get<0>(specE[Z]) += E;
	    std::get<1>(specE[Z]) += N;
	  }

	  if (COLL_SPECIES) {
	    speciesKey k1, k2;
	    CollCounts N;
	    while (1) {
	      MPI_Recv(&k1.first,  1, MPI_UNSIGNED_SHORT, i, 346, MPI_COMM_WORLD,
		       MPI_STATUS_IGNORE);
	      if (k1.first==255) break;
	      MPI_Recv(&k1.second, 1, MPI_UNSIGNED_SHORT, i, 347, MPI_COMM_WORLD,
		       MPI_STATUS_IGNORE);
	      MPI_Recv(&k2.first,  1, MPI_UNSIGNED_SHORT, i, 348, MPI_COMM_WORLD,
		       MPI_STATUS_IGNORE);
	      MPI_Recv(&k2.second, 1, MPI_UNSIGNED_SHORT, i, 349, MPI_COMM_WORLD,
		       MPI_STATUS_IGNORE);
	      MPI_Recv(&N[0],      1, MPI_UNSIGNED_LONG,  i, 350, MPI_COMM_WORLD,
			 MPI_STATUS_IGNORE);
	      MPI_Recv(&N[1],      1, MPI_UNSIGNED_LONG,  i, 351, MPI_COMM_WORLD,
			 MPI_STATUS_IGNORE);

	      dKey k(k1, k2);
	      if (collCount[0].find(k) == collCount[0].end()) 
		collCount[0][k] = N;
	      else {
		collCount[0][k][0] += N[0];
		collCount[0][k][1] += N[1];
	      }
	    }
	  }

	}

	mass  += val1;
	tempM += val2;
	consE += val3;
	consG += val4;
	totlE += val5;
	tempE += val6;
	elecE += val7;

      } // end: myid==0

    } // end: numprocs

    if (mass>0.0) {
      tempM /= mass;
      tempE /= mass;
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
    if (use_elec<0) Collide::printSpecies(spec, tempM);
    else printSpeciesElectrons(spec, tempM);
  } else if (aType == Weight) {	// Call the weighted printSpecies version
    printSpeciesElectrons(spec, tempM);
    printSpeciesColl();
  } else {			// Call the trace fraction version
    printSpeciesTrace();
  }

}

void CollideIon::printSpeciesColl()
{
  if (COLL_SPECIES) {

    unsigned long sum = 0;
    for (auto i : collCount[0]) {
      for (auto j : i.second) sum += j;
    }

    if (sum) {

      ostringstream sout;
      sout << outdir << runtag << ".DSMC_spc_log";
      ofstream mout(sout.str().c_str(), ios::app);

      // Print the header
      //
      mout << std::left
	   << std::setw(12) << "Time"      << std::setw(19) << tnow  << std::endl
	   << std::setw(12) << "Temp(ion)" << std::setw(18) << tempM << std::endl
	   << std::setw(12) << "Temp(elc)" << std::setw(18) << tempE << std::endl
	   << std::endl << std::right
	   << std::setw(20) << "<Sp 1|Sp 2> "
	   << std::setw(14) << "Scatter"
	   << std::setw(18) << "Frac scat"
	   << std::setw(14) << "Inelast"
	   << std::setw(18) << "Frac inel"
	   << std::endl << std::right
	   << std::setw(20) << "------------- " << std::right
	   << std::setw(10) << "---------"      << std::right
	   << std::setw(18) << "------------"   << std::right
	   << std::setw(10) << "-------"        << std::right
	   << std::setw(18) << "------------"   << std::endl;

      for (auto i : collCount[0]) {
      std::ostringstream sout;
      speciesKey k1 = i.first.first;
      speciesKey k2 = i.first.second;
      sout << "<" << std::setw(2) << k1.first 
	   << "," << std::setw(2) << k1.second 
	   << "|" << std::setw(2) << k2.first
	   << "," << std::setw(2) << k2.second << "> ";

      mout << std::setw(20) << std::right << sout.str() 
	   << std::setw(14) << i.second[0]
	   << std::setw(18) << static_cast<double>(i.second[0])/sum
	   << std::setw(14) << i.second[1]
	   << std::setw(18) << static_cast<double>(i.second[1])/sum
	   << std::endl;
      }
      mout << std::string(86, '-') << std::endl;
    }

    for (auto s : collCount) s.clear();
  }
}

void CollideIon::electronGather()
{
  bool IDBG = false;

  if ((aType==Direct or aType==Weight) && use_elec >= 0) {

    std::vector<double> eVel, iVel;

    // Interate through all cells
    //
    pHOT_iterator itree(*c0->Tree());
    
    ee.clear();

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

      if (NTC_DIST) {
	for (auto q : qv) {
	  NTC::NTCitem::vcMap v = ntcdb[itree.Cell()->mykey].CrsVel(q);
	  ee[q].push_back(v[electronKey]);
	}
      }
    }

    // Accumulate from threads
    //
    std::vector<double> loss, keE, keI, mom;
    unsigned Ovr=0, Acc=0, Tot=0;    for (int t=0; t<nthrds; t++) {
      loss.insert(loss.end(), velER[t].begin(), velER[t].end());
      velER[t].clear();
      Ovr += elecOvr[t];
      Acc += elecAcc[t];
      Tot += elecTot[t];
      elecOvr[t] = elecAcc[t] = elecTot[t] = 0;
    }

    if (ExactE and DebugE) {
      for (int t=0; t<nthrds; t++) {
	mom.insert(mom.end(), momD[t].begin(), momD[t].end());
      }
    }

    if (KE_DEBUG) {
      for (int t=0; t<nthrds; t++) {
	keE.insert(keE.end(), keER[t].begin(), keER[t].end());
	keER[t].clear();
	keI.insert(keI.end(), keIR[t].begin(), keIR[t].end());
	keIR[t].clear();
      }
    }

    if (NTC_DIST) {

      for (int n=1; n<numprocs; n++) {

	if (myid == n) {

	  int base = 326;

	  for (auto j : ee) {

	    double   val = j.first;
	    unsigned num = j.second.size();
	    
	    MPI_Send(&num, 1, MPI_UNSIGNED, 0, base+0, MPI_COMM_WORLD);
	    MPI_Send(&val, 1, MPI_DOUBLE,   0, base+1, MPI_COMM_WORLD);
	    MPI_Send(&j.second[0], num, MPI_DOUBLE, 0, base+2, MPI_COMM_WORLD);
	    base += 3;
	  }
	    
	  unsigned zero = 0;
	  MPI_Send(&zero, 1, MPI_UNSIGNED, 0, base, MPI_COMM_WORLD);
	  
	} // END: process send to root
	
	if (myid==0) {
	  
	  std::vector<double> v;
	  unsigned num;
	  double val;
	
	  int base = 326;

	  while (1) {
	      
	    MPI_Recv(&num, 1,    MPI_UNSIGNED, n, base+0, MPI_COMM_WORLD, 
		     MPI_STATUS_IGNORE);

	    if (num==0) break;

	    MPI_Recv(&val, 1,    MPI_DOUBLE,   n, base+1, MPI_COMM_WORLD, 
		     MPI_STATUS_IGNORE);
	    
	    v.resize(num);
	      
	    MPI_Recv(&v[0], num, MPI_DOUBLE,   n, base+2, MPI_COMM_WORLD,
		       MPI_STATUS_IGNORE);
	      
	    ee[val].insert( ee[val].end(), v.begin(), v.end() );
	    
	    base += 3;
	    
	  } // Loop over quantiles
	  
	} // Root receive loop

	MPI_Barrier(MPI_COMM_WORLD);

      } // Node loop

      if (myid==0) {
	for (auto &u : ee) {
	  eeHisto[u.first] = ahistoDPtr(new AsciiHisto<double>(u.second, 20, 0.01));
	}
      }

    } // END: NTC_DIST

    std::ofstream dbg;
    if (IDBG) {
      std::ostringstream sout;
      sout << runtag << ".eGather." << myid;
      dbg.open(sout.str().c_str(), ios::out | ios::app);
      sout.str(""); sout << "---- Step " << this_step
			 << " Time=" << tnow << " ";
      dbg << std::setw(70) << std::setfill('-') << left << sout.str()
	  << std::endl << std::setfill(' ');
    }

    (*barrier)("CollideIon::electronGather: BEFORE Send/Recv loop",
	       __FILE__, __LINE__);
    
    unsigned eNum;

    for (int i=1; i<numprocs; i++) {

      if (i == myid) {

	MPI_Send(&(eNum=eVel.size()), 1, MPI_UNSIGNED, 0, 435, MPI_COMM_WORLD);

	if (IDBG) dbg << std::setw(16) << "eVel.size() = " << std::setw(10) << eNum;

	if (eNum) MPI_Send(&eVel[0], eNum, MPI_DOUBLE, 0, 436, MPI_COMM_WORLD);
	if (eNum) MPI_Send(&iVel[0], eNum, MPI_DOUBLE, 0, 437, MPI_COMM_WORLD);

	if (IDBG) dbg << " ... eVel and iVel sent" << std::endl;

	MPI_Send(&(eNum=loss.size()), 1, MPI_UNSIGNED, 0, 438, MPI_COMM_WORLD);

	if (IDBG) dbg << std::setw(16) << "loss.size() = " << std::setw(10) << eNum;

	if (eNum) MPI_Send(&loss[0], eNum, MPI_DOUBLE, 0, 439, MPI_COMM_WORLD);
	
	if (IDBG) dbg << " ... loss sent" << std::endl;

	if (KE_DEBUG) {
	  MPI_Send(&(eNum=keE.size()), 1, MPI_UNSIGNED, 0, 440, MPI_COMM_WORLD);
	  if (IDBG) dbg << std::setw(16) << "keE.size() = " << std::setw(10) << eNum;

	  if (eNum) MPI_Send(&keE[0], eNum, MPI_DOUBLE, 0, 441, MPI_COMM_WORLD);
	  if (IDBG) dbg << " ... keE sent" << std::endl;

	  MPI_Send(&(eNum=keI.size()), 1, MPI_UNSIGNED, 0, 442, MPI_COMM_WORLD);
	  if (IDBG) dbg << std::setw(16) << "keI.size() = " << std::setw(10) << eNum;

	  if (eNum) MPI_Send(&keI[0], eNum, MPI_DOUBLE, 0, 443, MPI_COMM_WORLD);
	  if (IDBG) dbg << " ... keI sent" << std::endl;
	}

	if (ExactE and DebugE) {
	  MPI_Send(&(eNum=mom.size()), 1, MPI_UNSIGNED, 0, 444, MPI_COMM_WORLD);
	  if (IDBG) dbg << std::setw(16) << "mom.size() = " << std::setw(10) << eNum;

	  if (eNum) MPI_Send(&mom[0], eNum, MPI_DOUBLE, 0, 445, MPI_COMM_WORLD);
	  if (IDBG) dbg << " ... mom sent" << std::endl;
	  }
      }

				// Root receives from Node i
      if (0 == myid) {

	std::vector<double> vTmp;
	
	MPI_Recv(&eNum, 1, MPI_UNSIGNED, i, 435, MPI_COMM_WORLD,
		 MPI_STATUS_IGNORE);
	
	if (IDBG) dbg << "recvd from " << std::setw(4) << i
		      << std::setw(16) << " eVel.size() = "
		      << std::setw(10) << eNum;

	if (eNum) {
	  vTmp.resize(eNum);
	  
	  MPI_Recv(&vTmp[0], eNum, MPI_DOUBLE, i, 436, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  eVel.insert(eVel.begin(), vTmp.begin(), vTmp.end());

	  if (IDBG) dbg << " ... eVel recvd";

	  MPI_Recv(&vTmp[0], eNum, MPI_DOUBLE, i, 437, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  iVel.insert(iVel.begin(), vTmp.begin(), vTmp.end());

	  if (IDBG) dbg << ", iVel recvd" << std::endl;
	} else {
	  if (IDBG) dbg << std::endl;
	}

	MPI_Recv(&eNum, 1, MPI_UNSIGNED, i, 438, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if (IDBG) dbg << "recvd from " << std::setw(4) << i
		      << std::setw(16) << " loss.size() = "
		      << std::setw(10) << eNum;

	if (eNum) {
	  vTmp.resize(eNum);

	  MPI_Recv(&vTmp[0], eNum, MPI_DOUBLE, i, 439, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  loss.insert(loss.end(), vTmp.begin(), vTmp.end());

	  if (IDBG) dbg << " ... loss recvd" << std::endl;
	} else {
	  if (IDBG) dbg << std::endl;
	}

	if (KE_DEBUG) {
	  MPI_Recv(&eNum, 1, MPI_UNSIGNED, i, 440, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  if (IDBG) dbg << "recvd from " << std::setw(4) << i
			<< std::setw(16) << " keE.size() = " 
			<< std::setw(10) << eNum;

	  if (eNum) {
	    vTmp.resize(eNum);

	    MPI_Recv(&vTmp[0], eNum, MPI_DOUBLE, i, 441, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    keE.insert(keE.end(), vTmp.begin(), vTmp.end());

	    if (IDBG) dbg << " ... keE recvd" << std::endl;
	  } else {
	    if (IDBG) dbg << std::endl;
	  }

	  MPI_Recv(&eNum, 1, MPI_UNSIGNED, i, 442, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  if (IDBG) dbg << "recvd from " << std::setw(4) << i
			<< std::setw(16) << " keI.size() = "
			<< std::setw(10) << eNum;

	  if (eNum) {
	    vTmp.resize(eNum);

	    MPI_Recv(&vTmp[0], eNum, MPI_DOUBLE, i, 443, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    keI.insert(keI.end(), vTmp.begin(), vTmp.end());

	    if (IDBG) dbg << " ... keI recvd" << std::endl;
	  } else {
	    if (IDBG) dbg << std::endl;
	  }
	}

	if (ExactE and DebugE) {
	  MPI_Recv(&eNum, 1, MPI_UNSIGNED, i, 444, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  if (IDBG) dbg << "recvd from " << std::setw(4) << i
			<< std::setw(16) << " mom.size() = "
			<< std::setw(10) << eNum;

	  if (eNum) {
	    vTmp.resize(eNum);

	    MPI_Recv(&vTmp[0], eNum, MPI_DOUBLE, i, 445, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    mom.insert(mom.end(), vTmp.begin(), vTmp.end());

	    if (IDBG) dbg << " ... mom recvd" << std::endl;
	  } else {
	    if (IDBG) dbg << std::endl;
	  }
	}

      } // end: myid=0

    } // end: process loop

    (*barrier)("CollideIon::electronGather: AFTER Send/Recv loop",
	       __FILE__, __LINE__);

    MPI_Reduce(&Ovr, &Ovr_s, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Acc, &Acc_s, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&Tot, &Tot_s, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    
    (*barrier)("CollideIon::electronGather: AFTER REDUCE loop",
	       __FILE__, __LINE__);

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

      if (keE.size()) {
	keEH = ahistoDPtr(new AsciiHisto<double>(keE, 20, 0.01));
      }

      if (keI.size()) {
	keIH = ahistoDPtr(new AsciiHisto<double>(keI, 20, 0.01));
      }

      if (mom.size()) {
	momH = ahistoDPtr(new AsciiHisto<double>(mom, 20, 0.01));
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
  if (keEH.get()) {
    out << std::endl
	<< std::string(53, '-')  << std::endl
	<< "-----Relative electron energy gain/loss -------------" << std::endl
	<< std::string(53, '-')  << std::endl;
    (*keEH)(out);
  }
  if (keIH.get()) {
    out << std::endl
	<< std::string(53, '-')  << std::endl
	<< "-----Relative ion energy gain/loss ------------------" << std::endl
	<< std::string(53, '-')  << std::endl;
    (*keIH)(out);
  }

  if (momH.get()) {
    out << std::endl
	<< std::string(53, '-')  << std::endl
	<< "-----Electron momentum difference ratio -------------" << std::endl
	<< std::string(53, '-')  << std::endl;
    (*momH)(out);
  }


  out << std::string(53, '-') << std::endl
      << "-----Electron NTC diagnostics------------------------" << std::endl
      << std::string(53, '-') << std::endl << std::left
      << std::setw(14) << " Over"      << std::setw(16) << Ovr_s << std::endl 
      << std::setw(14) << " Accepted"  << std::setw(16) << Acc_s << std::endl
      << std::setw(14) << " Total"     << std::setw(16) << Tot_s << std::endl
      << std::fixed;

  if (Tot_s>0) 
    out << std::setw(14) << " Ratio"     << std::setw(16) << static_cast<double>(Acc_s)/Tot_s << std::endl 
	<< std::setw(14) << " Fail"      << std::setw(16) << static_cast<double>(Ovr_s)/Tot_s << std::endl;

  out << std::string(53, '-') << std::endl << std::right;

  out << std::endl;

  if (eeHisto.size() > 0) {

    for (auto j : eeHisto) {
      
      if (j.second.get()) {
	out << std::endl << std::string(53, '-') << std::endl
	    << std::left << std::fixed
	    << " Quantile: " << j.first << std::endl
	    << std::string(53, '-') << std::endl
	    << std::left << std::scientific;
	(*j.second)(out);
	out << std::endl;
      }
    }
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


void CollideIon::printSpeciesElectrons
(std::map<speciesKey, unsigned long>& spec, double temp)
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

    // Make species list
    //
    for (auto k : spec) specZ.insert(k.first.first);

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
      dout << std::setw(wid) << std::right << "Cons_E"
	   << std::setw(wid) << std::right << "Ions_E"
	   << std::setw(wid) << std::right << "Comb_E";
      if (use_elec>=0) {
	dout << std::setw(wid) << std::right << "Temp_E"
	     << std::setw(wid) << std::right << "Elec_E"
	     << std::setw(wid) << std::right << "Cons_G"
	     << std::setw(wid) << std::right << "Totl_E";
	for (auto Z : specZ) {
	  std::ostringstream sout1, sout2, sout3, sout4, sout5, sout6;
	  sout1 << "Eion(" << Z << ")";
	  sout2 << "Nion(" << Z << ")";
	  sout3 << "Tion(" << Z << ")";
	  sout4 << "Eelc(" << Z << ")";
	  sout5 << "Nelc(" << Z << ")";
	  sout6 << "Telc(" << Z << ")";
	  dout << std::setw(wid) << std::right << sout1.str()
	       << std::setw(wid) << std::right << sout2.str()
	       << std::setw(wid) << std::right << sout3.str()
	       << std::setw(wid) << std::right << sout4.str()
	       << std::setw(wid) << std::right << sout5.str()
	       << std::setw(wid) << std::right << sout6.str();
	}
      }
      dout << std::endl;
      
      dout << "# " 
	   << std::setw(wid) << std::right << "--------"
	   << std::setw(wid) << std::right << "--------";
      for (spCountMapItr it=spec.begin(); it != spec.end(); it++)
	dout << setw(wid) << std::right << "--------";
      dout << std::setw(wid) << std::right << "--------"
	   << std::setw(wid) << std::right << "--------"
	   << std::setw(wid) << std::right << "--------";
      if (use_elec>=0) {
	dout << std::setw(wid) << std::right << "--------"
	     << std::setw(wid) << std::right << "--------"
	     << std::setw(wid) << std::right << "--------"
	     << std::setw(wid) << std::right << "--------";
	for (size_t z=0; z<specZ.size(); z++)
	  dout << std::setw(wid) << std::right << "--------"
	       << std::setw(wid) << std::right << "--------"
	       << std::setw(wid) << std::right << "--------"
	       << std::setw(wid) << std::right << "--------"
	       << std::setw(wid) << std::right << "--------"
	       << std::setw(wid) << std::right << "--------";
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

  const double Tfac = 2.0*UserTreeDSMC::Eunit/3.0 * amu  /
    UserTreeDSMC::Munit/boltz;

  dout << std::setw(wid) << std::right << consE
       << std::setw(wid) << std::right << totlE
       << std::setw(wid) << std::right << totlE + consE;
  if (use_elec>=0)
    dout << std::setw(wid) << std::right << tempE
	 << std::setw(wid) << std::right << elecE
	 << std::setw(wid) << std::right << consG
	 << std::setw(wid) << std::right << elecE + totlE + consE + consG;
  for (auto Z : specZ) {
    if (specI.find(Z) != specI.end()) {
      double E = std::get<0>(specI[Z]);
      double N = std::get<1>(specI[Z]);
      dout << std::setw(wid) << std::right << E
	   << std::setw(wid) << std::right << N
	   << std::setw(wid) << std::right << E/N*Tfac;
    } else {
      dout << std::setw(wid) << std::right << 0.0
	   << std::setw(wid) << std::right << 0.0
	   << std::setw(wid) << std::right << 0.0;
    }
    if (specE.find(Z) != specE.end()) {
      double E = std::get<0>(specE[Z]);
      double N = std::get<1>(specE[Z]);
      dout << std::setw(wid) << std::right << E
	   << std::setw(wid) << std::right << N
	   << std::setw(wid) << std::right << E/N*Tfac;
    } else {
      dout << std::setw(wid) << std::right << 0.0
	   << std::setw(wid) << std::right << 0.0
	   << std::setw(wid) << std::right << 0.0;
    }
  }
  dout << std::endl;
}

void CollideIon::processConfig()
{
  // Parse test algorithm features
  //
  Configuration cfg;
  std::string config(config0);

  // Ensure that the original config is used, unless explicited edited
  // by the user
  //
  if (restart) config = runtag + ".CollideIon.config.json";

  try {

    if ( !boost::filesystem::exists(config) ) {
      if (myid==0) std::cout << "CollideIon: can't find config file <" 
			     << config << ">, using defaults" << std::endl;
    } else {
      cfg.load(config, "JSON");
    }
    
    if (!cfg.property_tree().count("_description")) {
      cfg.property_tree().put("_description", 
			      "This is a test config database for CollideIon");
      time_t t = time(0);   // get current time
      struct tm * now = localtime( & t );
      std::ostringstream sout;
      sout << (now->tm_year + 1900) << '-' 
	   << (now->tm_mon + 1) << '-'
	   <<  now->tm_mday;
      cfg.property_tree().put("_date", sout.str());
    }
    
    ExactE = 
      cfg.entry<bool>("ENERGY_ES", "Enable the explicit energy conservation algorithm", false);

    DebugE = 
      cfg.entry<bool>("ENERGY_ES_DBG", "Enable explicit energy conservation checking", true);

    AlgOrth = 
      cfg.entry<bool>("ENERGY_ORTHO", "Add energy in orthogonal direction", false);

    TRACE_ELEC =
      cfg.entry<bool>("TRACE_ELEC", "Add excess energy directly to the electrons", false);

    TRACE_REAPPLY =
      cfg.entry<bool>("TRACE_REAPPLY", "Immediately add the COM energy loss for splitting to the scattering interaction", false);

    COLL_SPECIES =
      cfg.entry<bool>("COLL_SPECIES", "Print collision count by species for debugging", false);

    SECONDARY_SCATTER =
      cfg.entry<bool>("SECONDARY_SCATTER", "Scatter electron with its donor ion", false);

    TRACE_FRAC =
      cfg.entry<double>("TRACE_FRAC", "Add this fraction to electrons and rest to ions", 1.0f);

    SAME_ELEC_SCAT = 
      cfg.entry<bool>("SAME_ELEC_SCAT", "Only scatter electrons with the same donor-ion mass", false);

    SAME_IONS_SCAT = 
      cfg.entry<bool>("SAME_IONS_SCAT", "Only scatter ions with the same mass", false);

    SAME_INTERACT = 
      cfg.entry<bool>("SAME_INTERACT", "Only perform interactions with equal-mass particles", false);

    SAME_TRACE_SUPP = 
      cfg.entry<bool>("SAME_TRACE_SUPP", "Distribute energy equally to trace species", false);

    NOCOOL_ELEC = 
      cfg.entry<bool>("NOCOOL_ELEC", "Suppress distribution of energy to electrons when using NOCOOL", false);

    NOSHARE_ELEC = 
      cfg.entry<bool>("NOSHARE_ELEC", "Suppress distribution of ionization energy between electrons", false);

    CLONE_ELEC = 
      cfg.entry<bool>("CLONE_ELEC", "Clone energy of ionizing electron to newly created free electron", false);

    frost_warning = 
      cfg.entry<bool>("frost_warning", "Warn if energy lost is smaller than available energy", false);

    DEBUG_SL =
      cfg.entry<bool>("DEBUG_SL", "Enable verbose interaction selection diagnostics", false);

    DEBUG_CR =
      cfg.entry<bool>("DEBUG_CR", "Enable printing of relative cross sections and probabilities for interaction selection", false);

    DEBUG_NQ =
      cfg.entry<bool>("DEBUG_NQ", "Printing of cross section debug info for unequal species only", false);

    NO_DOF =
      cfg.entry<bool>("NO_DOF", "Suppress adjustment of electron speed based on degrees of freedom", true);

    NO_VEL =
      cfg.entry<bool>("NO_VEL", "Suppress adjustment of electron speed for equipartition equilibrium", false);

    NO_ION_E =
      cfg.entry<bool>("NO_ION_E", "Suppress energy loss from ionization", false);

    NO_FF_E =
      cfg.entry<bool>("NO_FF_E", "Suppress energy loss from free-free", false);

    KE_DEBUG =
      cfg.entry<bool>("KE_DEBUG", "Check energy bookkeeping for weighted algorithm", true);

    DEBUG_THRESH =
      cfg.entry<double>("DEBUG_THRESH", "Threshold for reporting energy conservation bookkeeping", 1.0e-9);

    RECOMB_IP =
      cfg.entry<bool>("RECOMB_IP", "Electronic binding energy is lost in recombination", false);

    CROSS_DBG =
      cfg.entry<bool>("CROSS_DBG", "Enable verbose cross-section value diagnostics", false);

    EXCESS_DBG =
      cfg.entry<bool>("EXCESS_DBG", "Enable check for excess weight counter in trace algorithm", false);

    NTC_DIST =
      cfg.entry<bool>("NTC_DIST", "Enable NTC full distribution for electrons", false);

    FloorEv =
      cfg.entry<double>("FloorEv", "Minimum energy for Coulombic elastic scattering cross section", 0.05f);

    minCollFrac =
      cfg.entry<double>("minCollFrac", "Minimum relative fraction for collisional excitation", -1.0f);

    Collide::numSanityStop =
      cfg.entry<bool>("collStop", "Stop simulation if collisions per step are over threshold", false);

    Collide::numSanityMax =
      cfg.entry<unsigned>("maxStop", "Threshold for simulation stop", 100000000u);

    Collide::numSanityMsg =
      cfg.entry<bool>("collMsg", "Report collisions over threshold value", false);

    Collide::numSanityVal =
      cfg.entry<unsigned>("collMin", "Minimum threshold for reporting", 10000000u);

    Collide::numSanityFreq =
      cfg.entry<unsigned>("collFreq", "Stride for collision reporting", 2000000u);

    // Update atomic weight databases IF ElctronMass is specified
    // using direct call to underlying boost::property_tree
    //
    if (cfg.property_tree().find("ElectronMass.value") !=
	cfg.property_tree().not_found())
      {
	double mass = cfg.property_tree().get<double>("ElectronMass.value");
	UserTreeDSMC::atomic_weights[0] = atomic_weights[0] = mass;
      }

    if (myid==0) {
      // cfg.display();
      cfg.save(runtag +".CollideIon.config", "JSON");
    }
  }
  catch (...) {
    if (myid==0) std::cerr << "Error parsing CollideIon config info"
			   << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 54);
  }
}
