/*
  Compute initial conditions for a uniform density gas

  o ICs for the "direct" DSMC algorithm need constant particle number
    super particles for each species with the atomic numbers written to
    the species definition file.

  o ICs for the "weight" DSMC algorithm need variable particle number
    super particles for each species with the atomic numbers and
    weight factors written to the species definition file.

  o ICs for the "trace" DSMC algorithm should have constant
    super-particles masses with mass fractions for each species
    written to the species definition file altong with the position
    index for each trace fraction in the dattrib vector.
*/

#include <iostream>
#include <iomanip>
#include <numeric>
#include <cassert>
#include <memory>
#include <limits>
#include <random>
#include <tuple>

#include <yaml-cpp/yaml.h>

#include "Particle.H"
#include "globalInit.H"
#include "Species.H"
#include "atomic_constants.H"
#include <cxxopts.H>

#include <errno.h>
#include <sys/stat.h>

// Random types
//
typedef std::shared_ptr<std::mt19937> gen_ptr;
typedef std::shared_ptr<std::uniform_real_distribution<> > uniform_ptr;
typedef std::shared_ptr<std::normal_distribution<> > normal_ptr;

// EXP library support
//
#include <localmpi.H>
#include <libvars.H>

// Reference to n-body globals
//
using namespace __EXP__;	

#ifdef DEBUG
#include <fenv.h>
#include <fpetrap.h>

//===========================================
// Handlers defined in exputil/stack.cc
//===========================================

extern void mpi_print_trace(const string& routine, const string& msg,
			    const char *file, int line);

extern void mpi_gdb_print_trace(int sig);

extern void mpi_gdb_wait_trace(int sig);

//===========================================
// A signal handler to trap invalid FP only
//===========================================

void set_fpu_invalid_handler(void)
{
  // Flag invalid FP results only, such as 0/0 or infinity - infinity
  // or sqrt(-1).
  //
  feenableexcept(FE_INVALID);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_print_trace);
}

//===========================================
// A signal handler to produce a traceback
//===========================================

void set_fpu_trace_handler(void)
{
  // Flag all FP errors except inexact
  //
  // fedisableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Flag invalid FP results only, such as 0/0 or infinity - infinity
  // or sqrt(-1).
  //
  feenableexcept(FE_INVALID);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_print_trace);
}

//===========================================
// A signal handler to produce stop and wait
//===========================================

void set_fpu_gdb_handler(void)
{
  // Flag all FP errors except inexact
  //
  // fedisableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // Flag invalid FP results only, such as 0/0 or infinity - infinity
  // or sqrt(-1).
  //
  feenableexcept(FE_INVALID);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
  signal(SIGFPE, mpi_gdb_wait_trace);
}

#endif

//===========================================
// Write ChiantiPy script
//===========================================

bool file_exists(const std::string& fileName)
{
  std::ifstream infile(fileName);
  return infile.good();
}

void writeScript(void)
{
  const char *py =
#include "genIonization.h"
    ;

  const std::string file("genIonization");

  if (not file_exists(file)) {
    std::ofstream out(file);
    out << py;
    if (chmod(file.c_str(), S_IWUSR | S_IRUSR | S_IXUSR | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH)) {
      perror("Error in chmod:");
    }
  }
}

double Lunit;
double Tunit;
double Vunit;
double Munit;

std::vector<double> LL;

PeriodicTable PT;

//
// std RNGs
//
gen_ptr     gen;
uniform_ptr uniform;
normal_ptr  normal;

// ION collide types
//
enum Itype { Hybrid, Trace, Weight, Direct };
std::map<std::string, Itype> Types
{ {"Hybrid", Hybrid}, {"Trace", Trace}, {"Weight", Weight}, {"Direct", Direct} };

// Use CHIANTI or ION for ionization-recombination equilibrium
//
bool use_chianti   = true;
bool use_init_file = false;
bool use_yaml      = true;


// Model types
//
enum Mtype { Uniform, Interface };
std::map<std::string, Mtype> Models
{ {"Uniform", Uniform}, {"Interface", Interface} };


/**
   Make Uniform temperature box of gas
*/
void InitializeUniform(std::vector<Particle>& p, std::vector<double>& mass, double molW, double Ecut,
                       std::vector< std::map<unsigned char, double> >& T, std::vector<double> &L,
		       Itype type, int sp, int ne, int ni, int nd)
{
  unsigned npart = p.size();
  double   rho   = mass[0]/(L[0]*L[1]*L[2]);
  
  std::cout << std::string(70, '-')                 << std::endl;
  if (type == Trace) {
    std::cout << "MolW:        "  << molW                << std::endl;
    std::cout << "T(ion):      "  << T[0][0]  << " K "   << std::endl;
    std::cout << "T(elec):     "  << T[0][1]  << " K "   << std::endl;
  } else if (T.size()>1) {
    for (auto v : T[0]) {
      std::ostringstream sout;
      sout << "Temp " << PT[v.first]->name() << ":";
      std::cout << std::left << std::setw(13) << sout.str() << v.second
		<< " K "   << std::endl;
    }
  } else {
    std::cout << "Temp:        "  << T[0][0]  << " K "   << std::endl;
  }
  std::cout << "Number:      "  << npart            << std::endl;
  std::cout << "Length unit: "  << Lunit << " cm"   << std::endl;
  std::cout << "Time unit:   "  << Tunit << " s"    << std::endl;
  std::cout << "Vel unit:    "  << Vunit << " cm/s" << std::endl;
  std::cout << "Mass unit:   "  << Munit << " g"    << std::endl;
  std::cout << std::string(70, '-')                 << std::endl;
  
  /*
    We want every d.o.f. to have 1/2 k_B T of energy.  For classic
    DSMC, every superparticle has m/mu(Z) real particles [in physical
    units where mu(Z) = atomic_masses[Z]*amu].  So the velocity factor
    is given by the equality: 

    3/2*N*k_B*T = 3/2*m*k_B*T/mu(Z) = 3/2*m*v^2 

    where N is the number of particles: N=m/mu(Z), or 

    v^2 = k_B*T/mu(Z)

    For electrons:
    
    3/2*N_e*k_B*T = 3/2*m*eta*k_B*T/mu(Z) = 3/2*m/mu(Z)*eta*m_e*v_e^2 
    
    v_e^2 = k_B*T/m_e
  */

  std::map<unsigned char, double> varI, varE;
  for (auto v : T[0]) {
    unsigned char Z = v.first;
    if (Z>0) {			// All except TRACE
      varI[Z] = sqrt((boltz*T[0][Z])/(PT[Z]->weight()*amu)) / Vunit; // Ion
      varE[Z] = sqrt((boltz*T[0][Z])/(PT[0]->weight()*amu)) / Vunit; // Electron
    } else {			// TRACE
      varI[Z] = sqrt((boltz*T[0][0])/(molW*amu))            / Vunit; // Fiducial particle
      varE[Z] = sqrt((boltz*T[0][1])/(PT[0]->weight()*amu)) / Vunit; // Electrons
    }
  }
  
  double ttemp = T[0][0];
  
  for (unsigned i=0; i<npart; i++) {
    
    double KE = 0.0;
    
    p[i].iattrib.resize(ni, 0);
    p[i].dattrib.resize(nd, 0);

    if (type == Trace) {
      
      for (unsigned k=0; k<3; k++) {
	p[i].pos[k] = L[k]*((*uniform)(*gen));
	p[i].vel[k] = varI[0]*((*normal)(*gen));
      }
      
      if (ne>=0) {
	double totE;
	do {
	  totE = 0.0;
	  for (int l=0; l<3; l++) {
	    double v = (*normal)(*gen);
	    totE += v*v;
	    p[i].dattrib[ne+l] = varE[0] * v;
	  }
	} while (totE > Ecut);
      }

    } else {
      
      KeyConvert kc(p[i].iattrib[0]);
      speciesKey sKey  = kc.getKey();
      unsigned short Z = sKey.first;
      // unsigned short C = sKey.second;
      
      for (unsigned k=0; k<3; k++) {
	p[i].pos[k] = L[k]*(*uniform)(*gen);
	p[i].vel[k] = varI[Z] * (*normal)(*gen);
	KE += p[i].vel[k] * p[i].vel[k];
      }
      
      if (ne>=0) {
	for (int l=0; l<3; l++) {
	  p[i].dattrib[ne+l] = varE[Z] * (*normal)(*gen);
	}
      }
      
      ttemp = T[0][Z];
    }
    
    if (p[i].dattrib.size()>0) p[i].dattrib[0] = ttemp;
    if (p[i].dattrib.size()>1) p[i].dattrib[1] = rho;

    if (type == Direct) {
      if (p[i].dattrib.size()>2) p[i].dattrib[2] = KE;
    }
  }
}

/**
   Make split temperature-density box of gas
*/
void InitializeInterface(std::vector<Particle>& p,
			 std::vector<double>& mass, double molW,
			 std::vector< std::map<unsigned char, double> >& T,
			 std::vector<double> &L, Itype type,
			 int sp, int ne, int ni, int nd)
{
  unsigned npart = p.size();
  std::vector<double> rho
  { mass[0]/(0.5*L[0]*L[1]*L[2]), mass[1]/(0.5*L[0]*L[1]*L[2]) };
  
  std::cout << std::string(70, '-')                 << std::endl;

  for (size_t wh=0; wh<2; wh++) {

    std::cout << "Density = " << rho[wh] << std::endl;
    if (type == Trace) {
      std::cout << "     T(ion):      "  << T[wh][0]  << " K "   << std::endl
		<< "     T(elec):     "  << T[wh][1]  << " K "   << std::endl;
    } else if (T[wh].size()>1) {
      for (auto v : T[wh]) {
	std::ostringstream sout;
	sout << "     Temp " << PT[v.first]->name() << ":";
      std::cout << std::left << std::setw(13) << sout.str() << v.second
		<< " K "   << std::endl;
      }
    } else {
      std::cout << "     Temp:        "  << T[wh][0]  << " K "   << std::endl;
    }
  }

  std::cout << "Number:      "  << npart            << std::endl;
  std::cout << "Length unit: "  << Lunit << " cm"   << std::endl;
  std::cout << "Time unit:   "  << Tunit << " s"    << std::endl;
  std::cout << "Vel unit:    "  << Vunit << " cm/s" << std::endl;
  std::cout << "Mass unit:   "  << Munit << " g"    << std::endl;
  std::cout << std::string(70, '-')                 << std::endl;
  
  /*
    We want every d.o.f. to have 1/2 k_B T of energy.  For classic
    DSMC, every superparticle has m/mu(Z) real particles [in physical
    units where mu(Z) = atomic_masses[Z]*amu].  So the velocity factor
    is given by the equality: 

    3/2*N*k_B*T = 3/2*m*k_B*T/mu(Z) = 3/2*m*v^2 

    where N is the number of particles, or 

    v^2 = k_B*T/mu(Z)
  */
  std::vector< std::map<unsigned char, double> > varI(2), varE(2);
  for (unsigned char wh=0; wh<2; wh++) {
    for (auto v : T[wh]) {
      if (type == Trace) {
	if (v.first==0)		// Ion
	  varI[wh][0] = sqrt((boltz*v.second)/(molW*amu))            / Vunit;
	else			// Electron
	  varE[wh][0] = sqrt((boltz*v.second)/(PT[0]->weight()*amu)) / Vunit;
      } else {
	unsigned char Z = v.first;
	varI[wh][Z] = sqrt((boltz*T[wh][Z])/(PT[Z]->weight()*amu)) / Vunit; // Ion
	varE[wh][Z] = sqrt((boltz*T[wh][Z])/(PT[0]->weight()*amu)) / Vunit; // Electron
      } 
    }
  }
  
  double ttemp;

  std::uniform_int_distribution<> dist(0, 1);

  for (unsigned i=0; i<npart; i++) {
    
    double KE = 0.0;
    
    unsigned char wh = static_cast<unsigned char>(dist(*gen));


    p[i].iattrib.resize(ni, 0);
    p[i].dattrib.resize(nd, 0);

    if (type == Trace) {
      
      for (unsigned k=0; k<3; k++) {
	if (k==0)
	  p[i].pos[k] = 0.5*L[k]*((*uniform)(*gen) + wh);
	else
	  p[i].pos[k] = L[k]*(*uniform)(*gen);
	p[i].vel[k] = varI[wh][0]*(*normal)(*gen);
      }
      
      if (ne>=0) {
	for (int l=0; l<3; l++) {
	  p[i].dattrib[ne+l] = varE[wh][0] * (*normal)(*gen);
	}
      }

      ttemp = T[wh][0];

    } else {
      
      KeyConvert kc(p[i].iattrib[0]);
      speciesKey sKey  = kc.getKey();
      unsigned short Z = sKey.first;
      // unsigned short C = sKey.second;
      
      for (unsigned k=0; k<3; k++) {
	if (k==0)
	  p[i].pos[k] = 0.5*L[k]*((*uniform)(*gen) + wh);
	else
	  p[i].pos[k] = L[k]*(*uniform)(*gen);
	p[i].vel[k] = varI[wh][Z] * (*normal)(*gen);
	KE += p[i].vel[k] * p[i].vel[k];
      }
      
      if (ne>=0) {
	for (int l=0; l<3; l++) {
	  p[i].dattrib[ne+l] = varE[wh][Z] * (*normal)(*gen);
	}
      }
      
      ttemp = T[wh][Z];
    }
    
    if (p[i].dattrib.size()>0) p[i].dattrib[0] = ttemp;
    if (p[i].dattrib.size()>1) p[i].dattrib[1] = rho[wh];
    if (p[i].dattrib.size()>2) p[i].dattrib[2] = KE;
  }

}

void writeParticles(std::vector<Particle>& particles, const string& file, Itype type, 
		    const std::vector<double>& sF, 
		    const std::vector< std::vector<double> >& sI)
{
  // For tabulating mass fractions . . . 
  typedef std::tuple<double, unsigned> Tspc;
  typedef std::map<speciesKey, Tspc>   Frac;
  Frac frac;
  
  std::ofstream out(file.c_str());
  
  out.precision(10);
  
  out << setw(15) << particles.size()
      << setw(10) << particles[0].iattrib.size()
      << setw(10) << particles[0].dattrib.size()
      << std::endl;
  
  for (unsigned n=0; n<particles.size(); n++) {
    
    out << setw(18) << particles[n].mass;
    
    for (int k=0; k<3; k++) 
      out << setw(18) << particles[n].pos[k];
    
    for (int k=0; k<3; k++) 
      out << setw(18) << particles[n].vel[k];
    
    if (type != Trace) {
      KeyConvert kc(particles[n].iattrib[0]);
      speciesKey k = kc.getKey();
      Frac::iterator it = frac.find(k);
      if (it==frac.end()) frac[k] = Tspc(0.0, 0);
      std::get<0>(frac[k]) += particles[n].mass;
      std::get<1>(frac[k]) ++;
    }
    for (unsigned k=0; k<particles[n].iattrib.size(); k++)
      out << setw(12) << particles[n].iattrib[k];
    
    for (unsigned k=0; k<particles[n].dattrib.size(); k++)
      out << setw(18) << particles[n].dattrib[k];
    
    out << std::endl;
  }
  
  if (type != Trace) {
    
    double Mtot = 0.0;
    for (auto i : frac) Mtot += std::get<0>(i.second);
    
    double vol = 1.0; for (auto v : LL) vol *= v*Lunit;

    if ( type != Hybrid) {
      std::cout << std::setw( 3) << "Z"
		<< std::setw( 3) << "C"
		<< std::setw(16) << "Mass"
		<< std::setw(16) << "Fraction"
		<< std::setw(16) << "Target"
		<< std::setw(12) << "Count"
		<< std::endl
		<< std::setw( 3) << "-"
		<< std::setw( 3) << "-"
		<< std::setw(16) << "--------"
		<< std::setw(16) << "--------"
		<< std::setw(16) << "--------"
		<< std::setw(12) << "--------"
		<< std::endl;
      for (auto i : frac)
	std::cout << std::setw( 3) << i.first.first
		  << std::setw( 3) << i.first.second
		  << std::setw(16) << std::get<0>(i.second)
		  << std::setw(16) << std::get<0>(i.second)/Mtot
		  << std::setw(16) << sF[i.first.first-1] * sI[i.first.first-1][i.first.second-1]
		  << std::setw(12) << std::get<1>(i.second)
		  << std::endl;
    } else {
      std::cout << std::setw( 3) << "Z"
		<< std::setw( 3) << "C"
		<< std::setw(16) << "Mass"
		<< std::setw(16) << "Fraction"
		<< std::setw(12) << "Count"
		<< std::setw(16) << "n_Z (#/cc)"
		<< std::endl
		<< std::setw( 3) << "-"
		<< std::setw( 3) << "-"
		<< std::setw(16) << "--------"
		<< std::setw(16) << "--------"
		<< std::setw(12) << "--------"
		<< std::setw(16) << "----------"
		<< std::endl;
      for (auto i : frac)
	std::cout << std::setw( 3) << i.first.first
		  << std::setw( 3) << i.first.second
		  << std::setw(16) << std::get<0>(i.second)
		  << std::setw(16) << std::get<0>(i.second)/Mtot
		  << std::setw(12) << std::get<1>(i.second)
		  << std::setw(16) << std::get<0>(i.second)*Munit/(PT[i.first.first]->weight()*mp*vol)
		  << std::endl;
    }
    
    std::cout << std::string(70, '-') << std::endl
	      << "Empirical density (amu/cc) = " << Mtot*Munit/(mp*vol)
	      << std::endl << std::string(70, '-') << std::endl;
  }
}

void InitializeSpeciesDirect
(std::vector<Particle> & particles, 
 std::vector<unsigned char>& sZ, 
 std::vector<double>& sF, 
 std::vector< std::vector<double> >& sI,
 std::vector<double>& M,
 std::vector< std::map<unsigned char, double> >& T,
 int sp, int ne)
{
  std::vector< std::vector<double> > frac, cuml;
  
  //
  // Generate the ionization-fraction input file
  //
  for (auto n : sZ) {
    
    if (use_chianti) {
      
      const std::string ioneq("makeIonIC.ioneq");
      std::ostringstream sout;
      sout << "./genIonization"
	   << " -1 " << static_cast<unsigned>(n)
	   << " -2 " << static_cast<unsigned>(n)
	   << " -T " << T[0][n] << " -o " << ioneq;
      
      int ret = system(sout.str().c_str());
      
      if (ret) {
	std::cout << "System command  = " << sout.str() << std::endl;
	std::cout << "System ret code = " << ret << std::endl;
      }
      
      typedef std::vector<std::string> vString;
      
      std::string   inLine;
      std::ifstream sFile(ioneq.c_str());
      
      if (sFile.is_open()) {
	
	std::getline(sFile, inLine); // Read and discard the headers
	std::getline(sFile, inLine);
	
	{
	  vString s;
	  std::getline(sFile, inLine);
	  
	  std::istringstream iss(inLine);
	  std::copy(std::istream_iterator<std::string>(iss), 
		    std::istream_iterator<std::string>(), 
		    std::back_inserter<vString>(s));
	  
	  std::vector<double> v;
	  for (auto i : s) v.push_back(::atof(i.c_str()));
	  frac.push_back(v);
	}
	
	{
	  vString s;
	  std::getline(sFile, inLine);
	  
	  std::istringstream iss(inLine);
	  std::copy(std::istream_iterator<std::string>(iss), 
		    std::istream_iterator<std::string>(), 
		    std::back_inserter<vString>(s));
	  
	  std::vector<double> v;
	  for (vString::iterator i=s.begin(); i!=s.end(); i++)
	    v.push_back(::atof(i->c_str()));
	  cuml.push_back(v);
	}
      }
    }
    else {
      
      std::ostringstream sout;
      sout << "mpirun -np 1 genIonRecomb"
	   << " -Z " << static_cast<unsigned>(n)
	   << " -T " << T[0][n];
      
      int ret = system(sout.str().c_str());
      
      if (ret) {
	std::cout << "System command  = " << sout.str() << std::endl;
	std::cout << "System ret code = " << ret << std::endl;
      }
      
      typedef std::vector<std::string> vString;
      
      std::string   inLine;
      std::string   fName("IonRecombFrac.data");
      std::ifstream sFile(fName);
      
      if (sFile.is_open()) {
	
	vString s;
	std::getline(sFile, inLine);
	
	std::istringstream iss(inLine);
	std::copy(std::istream_iterator<std::string>(iss), 
		  std::istream_iterator<std::string>(), 
		  std::back_inserter<vString>(s));
	
	std::vector<double> v;
	for (auto i : s) v.push_back(::atof(i.c_str()));
	frac.push_back(v);
	
	double norm = 0.0;
	for (auto i : v) norm += i;
	
	if (fabs(norm - 1.0) > 1.0e-10) {
	  std::cout << "Normalization error: ";
	  for (auto i : v) std::cout << std::setw(16) << i;
	  std::cout << std::endl;
	}
	
	std::vector<double> c = v;
	for (size_t i=1; i<c.size(); i++) c[i] += c[i-1];
	cuml.push_back(c);
      } else {
	std::cerr << "Error opening ionization/recombination file "
		  << "<" << fName << "> for Z=" << n << std::endl;
	exit(-1);
      }
    }
  }
  
  int N = particles.size();
  
  // Compute cumulative species distribution
  //
  size_t NS = sF.size();

  // Normalize sF
  //
  double norm = std::accumulate(sF.begin(), sF.end(), 0.0);
  if (fabs(norm - 1.0)>1.0e-16) {
    std::cout << "Normalization change: " << norm << std::endl;
  }
  
  std::vector<double> frcS(NS), cumS(NS);
  for (size_t i=0; i<NS; i++) {
    sF[i]  /= norm;
    frcS[i] = sF[i]/PT[sZ[i]]->weight();
    cumS[i] = frcS[i] + (i ? cumS[i-1] : 0);
  }
  
  double normC = cumS.back();
  
  for (size_t i=0; i<NS; i++) {
    frcS[i] /= normC;
    cumS[i] /= normC;
  }
  
  double tKEi  = 0.0;
  double tKEe  = 0.0;
  double numbI = 0.0;
  double numbE = 0.0;
  double Eunit = Munit*Vunit*Vunit;

  for (int i=0; i<N; i++) {
    
    double rz = (*uniform)(*gen);
    double rc = (*uniform)(*gen);
    
    unsigned char Ci = 1, Zi;
    size_t indx;

    // Get the species
    for (indx=0; indx<NS; indx++) { 
      if (rz < cumS[indx]) {
	Zi = sZ[indx]; 
	break;
      }
    }
    // Get the ionization state
    for (int j=0; j<Zi+1; j++) {
      if (rc < cuml[indx][j]) {
	Ci = j+1;
	break;
      }
    }
    
    particles[i].mass  = M[0]/N * PT[sZ[indx]]->weight() * normC;
    
    KeyConvert kc(speciesKey(Zi, Ci));
    particles[i].iattrib[0] = kc.getInt();

    double KEi = 0.0, KEe = 0.0;
    
    for (int k=0; k<3; k++) {
      KEi += particles[i].vel[k] * particles[i].vel[k];
      if (ne>=0) {
	KEe += particles[i].dattrib[ne+k] * particles[i].dattrib[ne+k];
      }
    }
      
    if (particles[i].dattrib.size()>2) particles[i].dattrib[2] = KEi;

    // Ion KE
    //
    tKEi += 0.5 * particles[i].mass * KEi * Eunit;

    // Electron KE
    //
    tKEe += 0.5 * particles[i].mass * PT[0]->weight()/PT[Zi]->weight() * KEe * Eunit;

    // Ion number
    //
    numbI += particles[i].mass/PT[Zi]->weight() * Munit / amu;

    // Electron number
    //
    numbE += particles[i].mass/PT[Zi]->weight() * Munit / amu * Ci;
  }
  
  std::cout << "T (ion):     " << tKEi/(1.5*numbI*boltz) << std::endl
	    << "T (elec):    " << tKEe/(1.5*numbI*boltz) << std::endl
	    << std::string(70, '-') << std::endl;

  if (use_yaml) {
    YAML::Emitter out;

    out << YAML::BeginMap
	<< YAML::Key   << "species_map"
	<< YAML::BeginMap
	<< YAML::Key   << "type"
	<< YAML::Value << "direct"
	<< YAML::Key   << "elec"
	<< YAML::Value << sp
	<< YAML::Key   << "elements"
	<< YAML::Value << YAML::BeginSeq << YAML::Flow;
    for (auto v : sZ) out << v;
    out << YAML::EndSeq
	<< YAML::EndMap
	<< YAML::EndMap;
  
    std::ofstream fout("species.yml");

    fout << out.c_str() << std::endl;
    
  } else {
    std::ofstream out("species.spec");
    out << "direct" << std::endl;
    out << std::setw(6) << sp << std::endl;
    for (size_t indx=0; indx<NS; indx++) { 
      out << std::setw(6) << static_cast<unsigned>(sZ[indx])
	  << std::endl;
    }
  }
  
  sI = cuml;
  for (auto &s : sI) {
    double l = 0.0;
    for (auto &c : s) {
      double t = c;
      c -= l;
      l  = t;
    }
  }
}

void InitializeSpeciesWeight
(std::vector<Particle> & particles,
 std::vector<unsigned char>& sZ,
 std::vector<double>& sF,
 std::vector< std::vector<double> >& sI,
 std::vector<double>& M,
 std::vector< std::map<unsigned char, double> >& T,
 int sp, int ne)
{
  std::vector< std::vector<double> > frac, cuml;
  
  //
  // Generate the ionization-fraction input file
  //
  for (auto n : sZ) {
    
    if (use_chianti) {
      
      const std::string ioneq("makeIonIC.ioneq");
      std::ostringstream sout;
      sout << "./genIonization"
	   << " -1 " << static_cast<unsigned>(n)
	   << " -2 " << static_cast<unsigned>(n)
	   << " -T " << T[0][n] << " -o n" << ioneq;
      
      int ret = system(sout.str().c_str());
      
      if (ret) {
	std::cout << "System command  = " << sout.str() << std::endl;
	std::cout << "System ret code = " << ret << std::endl;
      }
      
      typedef std::vector<std::string> vString;
      
      std::string   inLine;
      std::ifstream sFile(ioneq.c_str());
      
      if (sFile.is_open()) {
	
	std::getline(sFile, inLine); // Read and discard the headers
	std::getline(sFile, inLine);
	
	{
	  vString s;
	  std::getline(sFile, inLine);
	  
	  std::istringstream iss(inLine);
	  std::copy(std::istream_iterator<std::string>(iss), 
		    std::istream_iterator<std::string>(), 
		    std::back_inserter<vString>(s));
	  
	  std::vector<double> v;
	  for (auto i : s) v.push_back(::atof(i.c_str()));
	  frac.push_back(v);
	}
	
	{
	  vString s;
	  std::getline(sFile, inLine);
	  
	  std::istringstream iss(inLine);
	  std::copy(std::istream_iterator<std::string>(iss), 
		    std::istream_iterator<std::string>(), 
		    std::back_inserter<vString>(s));
	  
	  std::vector<double> v;
	  for (vString::iterator i=s.begin(); i!=s.end(); i++)
	    v.push_back(::atof(i->c_str()));
	  cuml.push_back(v);
	}
      }
    }
    else {
      
      std::ostringstream sout;
      sout << "mpirun -np 1 genIonRecomb"
	   << " -Z " << static_cast<unsigned>(n)
	   << " -T " << T[0][n];
      
      int ret = system(sout.str().c_str());
      
      if (ret) {
	std::cout << "System command  = " << sout.str() << std::endl;
	std::cout << "System ret code = " << ret << std::endl;
      }
      
      typedef std::vector<std::string> vString;
      
      std::string   inLine;
      std::string   fName("IonRecombFrac.data");
      std::ifstream sFile(fName);
      
      if (sFile.is_open()) {
	
	vString s;
	std::getline(sFile, inLine);
	
	std::istringstream iss(inLine);
	std::copy(std::istream_iterator<std::string>(iss), 
		  std::istream_iterator<std::string>(), 
		  std::back_inserter<vString>(s));
	
	std::vector<double> v;
	for (auto i : s) v.push_back(::atof(i.c_str()));
	frac.push_back(v);
	
	double norm = 0.0;
	for (auto i : v) norm += i;
	
	if (fabs(norm - 1.0) > 1.0e-10) {
	  std::cout << "Normalization error: ";
	  for (auto i : v) std::cout << std::setw(16) << i;
	  std::cout << std::endl;
	}
	
	std::vector<double> c = v;
	for (size_t i=1; i<c.size(); i++) c[i] += c[i-1];
	cuml.push_back(c);
      } else {
	std::cerr << "Error opening ionization/recombination file "
		  << "<" << fName << "> for Z=" << n << std::endl;
	exit(-1);
      }
    }
  }
  
  int N = particles.size();
  // Compute cumulative species
  // distribution
  size_t NS = sF.size();
  // Normalize sF
  //
  double norm = std::accumulate(sF.begin(), sF.end(), 0.0);
  if (fabs(norm - 1.0)>1.0e-16) {
    std::cout << "Normalization change: " << norm << std::endl;
  }
  
  std::vector<double> frcS(sF), wght(NS);
  double fH = sF[0], W_H = 1.0;
  for (auto &v : frcS) v /= fH;
  
  wght[0] = W_H;
  for (size_t i=1; i<NS; i++)
    wght[i] = frcS[i]/PT[sZ[i]]->weight();
  
  std::uniform_int_distribution<> dist(0, NS-1);
  
  for (int i=0; i<N; i++) {
    // Get the species
    size_t indx = dist(*gen);
    unsigned char Ci = 1, Zi = sZ[indx];
    double rc = (*uniform)(*gen);
    
    // Get the ionization state
    for (int j=0; j<Zi+1; j++) {
      if (rc < cuml[indx][j]) {
	Ci = j+1;
	break;
      }
    }
    
    particles[i].mass  = M[0]/N * sF[indx] * NS;
    
    KeyConvert kc(speciesKey(Zi, Ci));
    particles[i].iattrib[0] = kc.getInt();
  }
  
  if (use_yaml) {
    YAML::Emitter out;

    out << YAML::BeginMap << YAML::Key << "species_map" << YAML::BeginMap;
    out << YAML::Key   << "type";
    out << YAML::Value << "weight";
    out << YAML::Key   << "cons";
    out << YAML::Value << sp;
    out << YAML::Key   << "elec";
    out << YAML::Value << ne;
    out << YAML::Key   << "elements";
    out << YAML::Value << YAML::BeginSeq;
    for (size_t indx=0; indx<NS; indx++) { 
      out << YAML::Flow << YAML::BeginSeq 
	  << static_cast<unsigned>(sZ[indx])
	  << wght[indx]
	  << M[0]/N * sF[indx] * NS
	  << YAML::EndSeq;
    }
    out << YAML::EndSeq;
    out << YAML::EndMap;
    out << YAML::EndMap;
    
    std::ofstream fout("species.yml");

    std::cout << out.c_str() << std::endl;

  } else {
    std::ofstream out("species.spec");
  
    out << "weight" << std::endl;
    out << std::setw(6) << sp;
    if (ne>=0) out << std::setw(6) << ne;
    out << std::endl;
  
    for (size_t indx=0; indx<NS; indx++) { 
      out << std::setw(6)  << static_cast<unsigned>(sZ[indx])
	  << std::setw(16) << wght[indx]
	  << std::setw(16) << M[0]/N * sF[indx] * NS
	  << std::endl;
    }
  }
  
  sI = cuml;
  for (auto &s : sI) {
    double l = 0.0;
    for (auto &c : s) {
      double t = c;
      c -= l;
      l = t;
    }
  }
}

void InitializeSpeciesHybrid
(std::vector<Particle> & particles, 
 std::vector<unsigned char>& sZ, 
 std::vector<double>& sF, 
 std::map<unsigned char, unsigned char> sC,
 std::vector< std::vector<double> >& sI,
 std::vector<double>& M,
 std::vector< std::map<unsigned char, double> >& T,
 int sp, int ne)
{
  std::map< unsigned short, std::vector<double> > frac;
  std::vector< std::vector<double> > cuml;
  
  //
  // Generate the ionization-fraction input file
  //
  for (auto n : sZ) {
    
    if (use_chianti) {
      
      const std::string ioneq("makeIonIC.ioneq");
      std::ostringstream sout;
      sout << "./genIonization"
	   << " -1 " << static_cast<unsigned>(n)
	   << " -2 " << static_cast<unsigned>(n)
	   << " -T " << T[0][n] << " -o " << ioneq;
      
      int ret = system(sout.str().c_str());
      
      if (ret) {
	std::cout << "System command  = " << sout.str() << std::endl;
	std::cout << "System ret code = " << ret << std::endl;
      }
      
      typedef std::vector<std::string> vString;
      
      std::string   inLine;
      std::ifstream sFile(ioneq.c_str());
      std::getline(sFile, inLine); // Read and discard the initial header
      
      if (sFile.is_open()) {
	
	std::getline(sFile, inLine);
	// Get the atomic species
	unsigned short Z;
	{
	  std::istringstream iss(inLine);
	  iss >> Z;
	}
	// Get the ionization fractions
	{
	  vString s;
	  std::getline(sFile, inLine);
	  
	  std::istringstream iss(inLine);
	  std::copy(std::istream_iterator<std::string>(iss), 
		    std::istream_iterator<std::string>(), 
		    std::back_inserter<vString>(s));
	  
	  std::vector<double> v;
	  for (auto i : s) v.push_back(::atof(i.c_str()));

	  if (sC.find(n) != sC.end()) {
	    std::vector<double> vv;
	    for (int i=0; i<sC[n]; i++) vv.push_back(v[i]);
	    v = vv;
	  }

	  frac[Z] = v;
	}
	
	{
	  vString s;
	  std::getline(sFile, inLine);
	  
	  std::istringstream iss(inLine);
	  std::copy(std::istream_iterator<std::string>(iss), 
		    std::istream_iterator<std::string>(), 
		    std::back_inserter<vString>(s));
	  
	  std::vector<double> v;
	  for (vString::iterator i=s.begin(); i!=s.end(); i++)
	    v.push_back(::atof(i->c_str()));

	  if (sC.find(n) != sC.end()) {
	    std::vector<double> vv;
	    for (int i=0; i<sC[n]; i++) vv.push_back(v[i]);
	    v = vv;
	  }

	  cuml.push_back(v);
	}
      }
    }
    else {
      
      std::ostringstream sout;
      sout << "mpirun -np 1 genIonRecomb"
	   << " -Z " << static_cast<unsigned>(n)
	   << " -T " << T[0][n];
      
      int ret = system(sout.str().c_str());
      
      if (ret) {
	std::cout << "System command  = " << sout.str() << std::endl;
	std::cout << "System ret code = " << ret << std::endl;
      }
      
      typedef std::vector<std::string> vString;
      
      std::string   inLine;
      std::string   fName("IonRecombFrac.data");
      std::ifstream sFile(fName);
      
      if (sFile.is_open()) {
	
	vString s;
	std::getline(sFile, inLine);
	
	std::istringstream iss(inLine);
	std::copy(std::istream_iterator<std::string>(iss), 
		  std::istream_iterator<std::string>(), 
		  std::back_inserter<vString>(s));
	
	std::vector<double> v;
	for (auto i : s) v.push_back(::atof(i.c_str()));

	if (sC.find(n) != sC.end()) {
	  std::vector<double> vv;
	  for (int i=0; i<sC[n]; i++) vv.push_back(v[i]);
	  v = vv;
	}

	frac[n] = v;
	
	double norm = 0.0;
	for (auto i : v) norm += i;
	
	if (fabs(norm - 1.0) > 1.0e-10) {
	  std::cout << "Normalization error: ";
	  for (auto i : v) std::cout << std::setw(16) << i;
	  std::cout << std::endl;
	}
	
	std::vector<double> c = v;
	for (size_t i=1; i<c.size(); i++) c[i] += c[i-1];
	cuml.push_back(c);
      } else {
	std::cerr << "Error opening ionization/recombination file "
		  << "<" << fName << "> for Z=" << n << std::endl;
	exit(-1);
      }
    }
  }
  
  int N = particles.size();
  // Compute cumulative species
  // distribution
  size_t NS = sF.size();

  // Normalize sF
  //
  double norm = std::accumulate(sF.begin(), sF.end(), 0.0);
  if (fabs(norm - 1.0)>1.0e-16) {
    std::cout << "Normalization change: " << norm << std::endl;
  }
  
  std::vector<double> frcS(sF), wght(NS);
  double fH = sF[0], W_H = 1.0;
  for (auto & v : frcS) v /= fH;
  
  wght[0] = W_H;
  for (size_t i=1; i<NS; i++)
    wght[i] = frcS[i]/PT[sZ[i]]->weight();
  
  std::uniform_int_distribution<> dist(0, NS-1);
  
  for (int i=0; i<N; i++) {
    // Get the species
    size_t indx = dist(*gen);
    unsigned char Ci = 0, Zi = sZ[indx];

    particles[i].mass  = M[0]/N * sF[indx] * NS;
    
    KeyConvert kc(speciesKey(Zi, Ci));
    particles[i].iattrib[0] = kc.getInt();
  }
  
  if (use_yaml) {
    YAML::Emitter out;
    out << YAML::BeginMap
	<< YAML::Key     << "species_map"
	<< YAML::BeginMap
	<< YAML::Key     << "type"
	<< YAML::Value   << "hybrid"
	<< YAML::Key     << "cons"
	<< YAML::Value   << sp
	<< YAML::Key     << "spos"
	<< YAML::Value   << sp+1;
    if (ne>=0) {
      out << YAML::Key   << "Elec"
	  << YAML::Value << ne;
    }
    out << YAML::Key << "elements" << YAML::BeginSeq;
    for (size_t indx=0; indx<NS; indx++) { 
      out << YAML::Flow << YAML::BeginSeq 
	  << static_cast<unsigned>(sZ[indx])
	  << wght[indx]
	  << M[0]/N * sF[indx] * NS
	  << YAML::EndSeq;
    }

    out << YAML::EndSeq
	<< YAML::EndMap
	<< YAML::EndMap;

    std::ofstream fout("species.yml");
    fout << out.c_str() << std::endl;

  } else {
    std::ofstream out("species.spec");
  
    out << "hybrid" << std::endl;
    out << std::setw(6) << sp;
    out << std::setw(6) << sp+1;
    if (ne>=0) out << std::setw(6) << ne;
    out << std::endl;
    
    for (size_t indx=0; indx<NS; indx++) { 
      out << std::setw(6)  << static_cast<unsigned>(sZ[indx])
	  << std::setw(16) << wght[indx]
	  << std::setw(16) << M[0]/N * sF[indx] * NS
	  << std::endl;
    }
  }
  
  sI = cuml;
  for (auto &s : sI) {
    double l = 0.0;
    for (auto &c : s) {
      double t = c;
      c -= l;
      l = t;
    }
  }
}

void InitializeSpeciesTrace
(std::vector<Particle> & particles, 
 std::vector<unsigned char>& sZ, 
 std::vector<double>& sF,
 std::map<unsigned char, unsigned char> sC,
 std::vector<double>& M,
 std::vector< std::map<unsigned char, double> >& T,
 std::vector<double>& L,
 Mtype model, int sp, int ne, bool ECtype)
{
  size_t Ncomp = T.size();

  typedef std::vector< std::vector<double> > spA;
  std::vector<spA> frac(Ncomp), cuml(Ncomp);
  
  if (use_init_file and Ncomp==1) {
    std::ifstream sFile("IonRecombFrac.data");

    typedef std::vector<std::string> vString;
    std::string inLine;

    if (sFile.is_open()) {
      vString s;
      std::getline(sFile, inLine);
	
      std::istringstream iss(inLine);
      std::copy(std::istream_iterator<std::string>(iss), 
		std::istream_iterator<std::string>(), 
		std::back_inserter<vString>(s));
      
      std::vector<double> v;
      for (auto i : s) v.push_back(::atof(i.c_str()));
      
      int cnt = 0;
      for (auto n : sZ) {
	double sum = 0.0;
	std::vector<double> V;
	int Cmax = n;
	if (sC.find(n)!=sC.end()) Cmax = sC[n]-1;
	for (int C=0; C<Cmax; C++) {
	  sum += v[cnt];
	  V.push_back(v[cnt]);
	  cnt++;
	}
	V.push_back(1.0 - sum);
	frac[0].push_back(V);

	double norm = 0.0;
	for (auto i : V) norm += i;
	
	if (fabs(norm - 1.0) > 1.0e-10) {
	  std::cout << "Normalization error: ";
	  for (auto i : V) std::cout << std::setw(16) << i;
	  std::cout << std::endl;
	}
	
	std::vector<double> c = V;
	for (size_t i=1; i<c.size(); i++) c[i] += c[i-1];
	cuml[0].push_back(c);
      }
    }

  }
  else {

    //
    // Generate the ionization-fraction input file
    //
    for (size_t nc=0; nc<Ncomp; nc++) {

      for (auto n : sZ) {
      
	if (use_chianti) {
	
	  const std::string ioneq("makeIonIC.ioneq");
	  std::ostringstream sout;
	  sout << "./genIonization"
	    << " -1 " << static_cast<unsigned>(n)
	       << " -2 " << static_cast<unsigned>(n)
	       << " -T " << T[nc][1] << " -o " << ioneq;
	
	  int ret = system(sout.str().c_str());
	  
	  if (ret) {
	    std::cout << "System command  = " << sout.str() << std::endl;
	    std::cout << "System ret code = " << ret << std::endl;
	  }
	
	  typedef std::vector<std::string> vString;
	  
	  std::string   inLine;
	  std::ifstream sFile(ioneq.c_str());
	  if (sFile.is_open()) {
	    
	    std::getline(sFile, inLine); // Read and discard the headers
	    std::getline(sFile, inLine);
	  
	    {
	      vString s;
	      std::getline(sFile, inLine);
	      
	      std::istringstream iss(inLine);
	      std::copy(std::istream_iterator<std::string>(iss), 
			std::istream_iterator<std::string>(), 
			std::back_inserter<vString>(s));
	      
	      std::vector<double> v;
	      for (vString::iterator i=s.begin(); i!=s.end(); i++)
		v.push_back(::atof(i->c_str()));

	      if (sC.find(n) != sC.end()) {
		std::vector<double> vv;
		for (int i=0; i<sC[n]; i++) vv.push_back(v[i]);
		v = vv;
	      }

	      frac[nc].push_back(v);
	    }
	  
	    {
	      vString s;
	      std::getline(sFile, inLine);
	      
	      std::istringstream iss(inLine);
	      std::copy(std::istream_iterator<std::string>(iss), 
			std::istream_iterator<std::string>(), 
			std::back_inserter<vString>(s));
	      
	      std::vector<double> v;
	      for (vString::iterator i=s.begin(); i!=s.end(); i++)
		v.push_back(::atof(i->c_str()));

	      if (sC.find(n) != sC.end()) {
		std::vector<double> vv;
		for (int i=0; i<sC[n]; i++) vv.push_back(v[i]);
		v = vv;
	      }

	      cuml[nc].push_back(v);
	    }
	  }
	}
	else {
	  
	  std::ostringstream sout;
	  sout << "mpirun -np 1 genIonRecomb"
	       << " -Z " << static_cast<unsigned>(n)
	       << " -T " << T[nc][1];
	
	  int ret = system(sout.str().c_str());
      
	  if (ret) {
	    std::cout << "System command  = " << sout.str() << std::endl;
	    std::cout << "System ret code = " << ret << std::endl;
	  }
	
	  typedef std::vector<std::string> vString;
      
	  std::string   inLine;
	  std::ifstream sFile("IonRecombFrac.data");
	  if (sFile.is_open()) {
	    
	    vString s;
	    std::getline(sFile, inLine);
	    
	    std::istringstream iss(inLine);
	    std::copy(std::istream_iterator<std::string>(iss), 
		      std::istream_iterator<std::string>(), 
		      std::back_inserter<vString>(s));
	    
	    std::vector<double> v;
	    for (auto i : s) v.push_back(::atof(i.c_str()));

	    if (sC.find(n) != sC.end()) {
	      std::vector<double> vv;
	      for (int i=0; i<sC[n]; i++) vv.push_back(v[i]);
	      v = vv;
	    }

	    frac[nc].push_back(v);
	  
	    double norm = 0.0;
	    for (auto i : v) norm += i;
	    
	    if (fabs(norm - 1.0) > 1.0e-10) {
	      std::cout << "Normalization error: ";
	      for (auto i : v) std::cout << std::setw(16) << i;
	      std::cout << std::endl;
	    }
	  
	    std::vector<double> c = v;
	    for (size_t i=1; i<c.size(); i++) c[i] += c[i-1];
	    cuml[nc].push_back(c);
	  }
	}
      }
    } // END: number of components
  }
  
  int N = particles.size();
  
  // Compute cumulative species distribution
  //
  int NS = sF.size();

  double molW = 0.0;
  for (int i=0; i<NS; i++)
    molW += sF[i]/PT[sZ[i]]->weight();
  molW = 1.0/molW;
  
  // Normalize sF
  //
  double norm = std::accumulate(sF.begin(), sF.end(), 0.0);
  if (fabs(norm - 1.0)>1.0e-16) {
    std::cout << "Normalization change: " << norm << std::endl;
  }
  for (int i=0; i<NS; i++) sF[i]  /= norm;
  
  for (size_t nc=0; nc<Ncomp; nc++) {
    for (int indx=0; indx<NS; indx++) { 
      double norm = std::accumulate(frac[nc][indx].begin(), frac[nc][indx].end(), 0.0);
      for (auto &i : frac[nc][indx]) i /= norm;
    }
  }
  
  std::vector<double> eta(Ncomp, 1.0);

  double tKEi = 0.0, tKEe = 0.0, numb = 0.0;

  for (int i=0; i<N; i++) {
    
    size_t wh = 0;
    if (model == Interface) {
      if (particles[i].pos[0] > 0.5*L[0]) wh = 1;
    }

    particles[i].mass  = M[wh]/N;
    
    // Sanity check
    //
    double test = 0.0;

    int cur = sp;
    double Eta = 0.0, Mol = 0.0;
    // Get the species
    //
    for (int indx=0; indx<NS; indx++) { 
      // Get the ionization state
      //
      int maxC = sZ[indx]+1;	// Enforce upper ionization-state limit
      auto it = sC.find(sZ[indx]);
      if (it != sC.end()) maxC = it->second;

      double cc = sF[indx]/PT[sZ[indx]]->weight();
      Mol += cc;

      for (int j=0; j<maxC; j++) {
	particles[i].dattrib[cur++] = sF[indx]*frac[wh][indx][j];
	test += sF[indx] * frac[wh][indx][j];
	Eta  += frac[wh][indx][j] * cc * j;
      }
    }
    assert ( fabs(test-1.0) < 1.0e-12 );

    eta[wh] = Eta/Mol;

    double KEi = 0.0, KEe = 0.0;
    for (int k=0; k<3; k++) {
      KEi += particles[i].vel[k] * particles[i].vel[k];
      if (ne>=0) {
	if (ECtype)		// mean-mass correction
	  particles[i].dattrib[ne+k] *= 1.0/sqrt(eta[wh]);
	KEe += particles[i].dattrib[ne+k] * particles[i].dattrib[ne+k];
      }
    }
      
    if (particles[i].dattrib.size()>2) particles[i].dattrib[2] = KEi;
    
    // Kinetic energies
    //
    tKEi += 0.5*particles[i].mass * KEi;
    tKEe += 0.5*particles[i].mass * KEe * PT[0]->weight() * eta[wh]/ molW;

    // Ion number
    //
    numb += particles[i].mass/molW * Munit / amu;
  }
  
  if (use_yaml) {
    YAML::Emitter out;

    out << YAML::BeginMap
	<< YAML::Key << "species_map"
	<< YAML::BeginMap
	<< YAML::Key   << "type"
	<< YAML::Value << "trace"
	<< YAML::Key << "cons"
	<< YAML::Value << sp - 1
	<< YAML::Key << "elec"
	<< YAML::Value << ne
	<< YAML::Key << "elements"
	<< YAML::Value << YAML::BeginSeq;
  
    int cntr = sp;
    for (int indx=0; indx<NS; indx++) { 
      int maxC = sZ[indx]+1;	// Enforce upper ionization-state limit
      auto it = sC.find(sZ[indx]);
      if (it != sC.end()) maxC = it->second;

      for (int j=0; j<maxC; j++) {
	out << YAML::Flow << YAML::BeginSeq 
	    << static_cast<unsigned>(sZ[indx])
	    << j + 1
	    << cntr++ << YAML::EndSeq;
      }
    }
    
    out << YAML::EndSeq
	<< YAML::EndMap
	<< YAML::EndMap;

    std::ofstream fout("species.yml");
    fout << out.c_str() << std::endl;

  } else {

    std::ofstream out("species.spec");
    out << "trace" << std::endl;
    // Conservation position and electron position (-1 for none)
    //
    out << std::setw(6) << sp-1
	<< std::setw(6) << ne << std::endl;
    // Starting species position
    //
    int cntr = sp;
    for (int indx=0; indx<NS; indx++) { 
      int maxC = sZ[indx]+1;	// Enforce upper ionization-state limit
      auto it = sC.find(sZ[indx]);
      if (it != sC.end()) maxC = it->second;

      for (int j=0; j<maxC; j++) {
	out << std::setw(6) << static_cast<unsigned>(sZ[indx])
	    << std::setw(6) << j + 1
	    << std::setw(6) << cntr++
	    << std::endl;
      }
    }
  }
  
  double Eunit = Munit*Vunit*Vunit;

  std::cout << std::string(70, '-')    << std::endl
	    << "KE (ion):    " << tKEi << std::endl
	    << "KE (elec)    " << tKEe << std::endl
	    << std::string(70, '-')    << std::endl
	    << "MolW:        " << molW << std::endl
	    << "T (ion):     " << tKEi*Eunit/(1.5*numb*boltz) << std::endl
	    << "T (elec):    " << tKEe*Eunit/(1.5*numb*boltz) << std::endl
	    << std::string(70, '-')    << std::endl;

} // END: writeParticles

int
main (int ac, char **av)
{
  Itype type   = Direct;
  Mtype model  = Uniform;
  bool  ECtype = false;

  double D, L, R, Temp, Telec, Ecut;
  std::string config;
  std::string oname;
  unsigned seed;
  int ne = -1, ni = 2, nd = 5;
  int npart;
  
  std::string cmd_line;
  for (int i=0; i<ac; i++) {
    cmd_line += av[i];
    cmd_line += " ";
  }
  
  cxxopts::Options options(av[0], "Compute gas initial conditions for DSMC\n");

  options.add_options()
    ("h,help", "produce help message")
    ("electrons", "set up for trace, weighted or hybrid species with electrons")
    ("yaml", "write YAML species config file")
    ("old", "write old-style species config file")
    ("traceEC", "use explicit conservation interactions mode for Trace")
    ("C,CHIANTI", "use CHIANTI to set recombination-ionization equilibrium",
     cxxopts::value<bool>(use_chianti)->default_value("false"))
    ("I,INIT", "use init file to set recombination-ionization equilibrium",
     cxxopts::value<bool>(use_init_file)->default_value("false"))
    ("D,dens", "density in particles per cc",
     cxxopts::value<double>(D)->default_value("1.0"))
    ("T,temp", "override config file temperature for Trace, if >0",
     cxxopts::value<double>(Temp)->default_value("-1.0"))
    ("Telec", "temperature for electrons, if Telec>0",
     cxxopts::value<double>(Telec)->default_value("-1.0"))
    ("L,length", "box scale in system length units",
     cxxopts::value<double>(L)->default_value("1.0"))
    ("l,Lunit", "length in system units",
     cxxopts::value<double>(Lunit)->default_value("1.0"))
    ("t,Tunit", "time unit in years",
     cxxopts::value<double>(Tunit)->default_value("1.0e3"))
    ("m,Munit", "Mass unit in solar masses",
     cxxopts::value<double>(Munit)->default_value("1.0"))
    ("R,ratio", "slab length ratio (1 is cube)",
     cxxopts::value<double>(R)->default_value("1.0"))
    ("i,num-int", "number of integer attributes",
     cxxopts::value<int>(ni)->default_value("2"))
    ("d,num-double", "base number of double attributes",
     cxxopts::value<int>(nd)->default_value("6"))
    ("N,number", "number of particles",
     cxxopts::value<int>(npart)->default_value("10000"))
    ("E,Ecut", "truncation of electron tail in kT",
     cxxopts::value<double>(Ecut)->default_value("1.0e20"))
    ("c,config", "element config file",
     cxxopts::value<std::string>(config)->default_value("makeIon.config"))
    ("o,output", "output prefix",
     cxxopts::value<std::string>(oname)->default_value("out"))
    ;
  
  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    std::cout << "Example: " << std::endl;
    std::cout << "\t" << av[0]
	      << "--length=0.0001 --number=1000 --dens=0.01 --electrons -c trace.config --output=out"
	      << std::endl;
    return 1;
  }
  
  if (vm.count("yaml")) {
    use_yaml = true;
  }

  if (vm.count("old")) {
    use_yaml = false;
  }

  if (vm.count("traceEC")) {
    ECtype = true;
  }

  if (myid==0) {
    std::string prefix("makeIon");
    std::string cmdFile = prefix + "." + oname + ".cmd_line";
    std::ofstream out(cmdFile.c_str());
    if (!out) {
      std::cerr << "makeIon: error opening <" << cmdFile
                << "> for writing" << std::endl;
    } else {
      out << cmd_line << std::endl;
    }
  }
  
  if (use_chianti) writeScript();

#ifdef DEBUG                    // For gdb . . . 
  sleep(20);
  // set_fpu_handler();         // Make gdb trap FPU exceptions
  set_fpu_gdb_handler();	// Make gdb trap FPU exceptions
#endif

  Lunit  *= pc;
  Tunit  *= year;
  Munit  *= msun;
  Vunit   = Lunit/Tunit;
  
  gen     = std::make_shared<std::mt19937>(seed);
  uniform = std::make_shared<std::uniform_real_distribution<>>(0.0, 1.0);
  normal  = std::make_shared<std::normal_distribution<>>(0.0, 1.0);
  
  //
  // Define the atomic species statistics
  //
  
  // Element data
  std::vector<unsigned char> sZ;
  std::vector<double>        sF;

  // Maximum ionization level
  std::map<unsigned char, unsigned char> sC;

  // Temperatures
  std::vector< std::map<unsigned char, double> > T;
  std::vector<double> rho;

  // Species position
  int sp = -1;

  // Default molecular weight
  double molW = 1.0;

  // Parse element file
  {
    // Load the json/yaml file
    //
    YAML::Node iroot = YAML::LoadFile(config);
  
    double norm = 0.0;

    // Look for the type
    //
    std::string st("Direct");

    if (iroot["type"]) {
      st = iroot["type"].as<std::string>();
      
      if (Types.find(st) == Types.end()) {
	std::cout << "Type <" << st << "> is not valid" << std::endl
		  << "Valid types are:";
	for (auto v : Types) std::cout << " " << v.first;
	std::cout << std::endl;
	exit(-1);
      }
    } 
    
    type = Types[st];

    std::string md("Uniform");
    if (iroot["model"]) {
      md = iroot["model"].as<std::string>();

      if (Models.find(md) == Models.end()) {
	std::cout << "Model <" << md << "> is not valid" << std::endl
		  << "Valid types are:";
	for (auto v : Models) std::cout << " " << v.first;
	std::cout << std::endl;
	exit(-1);
      }
    }

    model = Models[md];

    if (iroot["electrons"]) {
      ne = 0;
    }
  
    if (model==Interface) {
      T.resize(2);
    } else {
      T.resize(1);
    }

    double fH = 0.0;

    if (iroot["elements"]) {

      for (YAML::const_iterator it=iroot["elements"].begin(); it != iroot["elements"].end(); ++it) {
	
	unsigned i = it->first.as<int>();

	sZ.push_back(i);

	for (YAML::const_iterator jt=it->second.begin(); jt != it->second.end(); ++jt) {
	  std::string    tag = jt->first.as<std::string>();
	  
	  if (tag == "logN") {
	    double val = jt->second.as<double>();
	    if (i==1) fH = val;
	    double wgt = pow(10.0, val - fH)*PT[i]->weight();
	    sF.push_back(wgt);
	    norm += sF.back();
	  } else if (tag == "mfrac") {
	    double val = jt->second.as<double>();
	    sF.push_back(val);
	    norm += sF.back();
	  } else if (tag == "cmax") {
	    unsigned char cval = jt->second.as<unsigned char>();
	    sC[i] = cval;
	  } else {
	    std::cerr << "Missing element definition for atomic number " << i
		      << std::endl;
	    exit(-3);
	  }
	}

	if (type != Trace) {
	  if (it->second["temp"]) {
	    T[0][i] = it->second["temp"].as<double>();
	  }
	}
      }
    }

    if (norm>0.0) {
      for (auto & v : sF) v /= norm;
    } else {
      std::cout << "Error: zero mass fraction norm" << std::endl;
      exit(-1);
    }

    if (type == Trace) {
      molW = 0.0;
      for (size_t i=0; i<sF.size(); i++)
	molW += sF[i]/PT[sZ[i]]->weight();
      molW = 1.0/molW;
    }

    if (model==Interface) {
      rho.resize(2);
      
      if (iroot["components"]) {
	
	// Iterate through component sequence
	for (YAML::const_iterator it=iroot["components"].begin(); it != iroot["elements"].end(); ++it) {
	
	  for (YAML::const_iterator jt=it->second.begin(); jt != it->second.end(); ++jt) {
	    unsigned j = jt->first.as<int>(); // Unpack component variables
	    if (j>0 and j<2) {
	      rho[j-1]  = jt->second["Density"].as<double>();
	      T[j-1][0] = jt->second["Temp"   ].as<double>();
	      T[j-1][1] = jt->second["Etemp"  ].as<double>();
	    }
	  }
	    
	}

      }
    } else {
      rho.push_back(D);
      if (Temp>0.0) T[0][0]  = Temp;
      else          T[0][0]  = iroot["temp"].as<double>();
      T[0][1] = T[0][0];
      if (Telec>0.0) T[0][1] = Telec;
      else           T[0][1] = iroot["telc"].as<double>();
    }
      
    if (type==Trace) {
      nd++;			// Conservation of energy position
      sp = nd;

				// Augment for species number
      for (size_t indx=0; indx<sF.size(); indx++) { 
	int maxC = sZ[indx]+1;	// Enforce upper ionization-state limit
	auto it = sC.find(sZ[indx]);
	if (it != sC.end()) maxC = it->second;

	for (int j=0; j<maxC; j++) nd++;
      }

      // Electron velocities and cons
      if (ne>=0) {
	ne = nd;		// Electron start position
	nd += 4;		// Electron velocity and conservation
      }

    }
    else if (type==Hybrid) {
      auto it = std::max_element(std::begin(sZ), std::end(sZ));
      size_t maxSp = *it;
      nd++;			// Energy conservation 
      sp = nd++;		// Species position
      nd += maxSp;
      if (ne>=0) {
	ne = nd;		// Electron start position
	nd += 4;		// Electron velocity and conservation
      }
    }
    else if (type==Weight or type==Direct) {
      sp = ++nd;		// Energy conservation 
      if (ne>=0) {
	ne = nd;		// Electron start position
	nd += 4;		// Electron velocity and conservation
      }
    }
  }

  // Ion fractions
  std::vector< std::vector<double> > sI;

  /* Additional species, e.g.
     sF[2] = 0.0; sZ[2] = 3;  // Li
     sF[3] = 0.0; sZ[3] = 6;  // C
     sF[4] = 0.0; sZ[4] = 7;  // N
     sF[5] = 0.0  sZ[5] = 8;  // O
     sF[6] = 0.0; sZ[6] = 12; // Mg
  */
  
  // Cube axes
  //
  LL.resize(3, L);
  LL[0] *= R;
  double vol = 1.0;
  for (auto v : LL) vol *= v*Lunit;
  
  // Mass in box in m_p
  //
  std::vector<double> Mass;
  
  std::vector<Particle> particles(npart);
  
  // Initialize the phase space vector
  //
  switch (model) {
  case Interface:
    Mass.push_back( mp*rho[0]*0.5*vol/Munit );
    Mass.push_back( mp*rho[1]*0.5*vol/Munit );
    InitializeInterface(particles, Mass, molW, T, LL, type, sp, ne, ni, nd);
    break;
  case Uniform:
  default:
    Mass.push_back(mp*D*vol/Munit);
    InitializeUniform(particles, Mass, molW, Ecut, T, LL, type, sp, ne, ni, nd);
    break;
  }

  // Initialize the Z, C's	
  //
  switch (type) {
  case Hybrid:
    InitializeSpeciesHybrid(particles, sZ, sF, sC, sI, Mass, T, sp, ne);
    break;
  case Weight:
    InitializeSpeciesWeight(particles, sZ, sF, sI, Mass, T, sp, ne);
    break;
  case Trace:
    InitializeSpeciesTrace (particles, sZ, sF, sC, Mass, T, LL, model, sp, ne,
			    ECtype);
    // Compute molecular weight
    molW = 0.0;
    for (size_t k=0; k<sZ.size(); k++) molW += sF[k]/PT[sZ[k]]->weight();
    molW = 1.0/molW;
    break;
  case Direct:
  default:
    InitializeSpeciesDirect(particles, sZ, sF, sI, Mass, T, sp, ne);
  }
  
  // Output file
  //
  writeParticles(particles, oname + ".bod", type, sF, sI);
  
  return 0;
}
