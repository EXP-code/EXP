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
#include <tuple>

#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <boost/shared_ptr.hpp>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "Particle.H"
#include "globalInit.H"
#include "Species.H"

//
// Boost random types
//
typedef boost::shared_ptr<boost::mt19937> gen_ptr;
typedef boost::shared_ptr<boost::uniform_real<> > uniform_ptr;
typedef boost::shared_ptr<boost::normal_distribution<> > normal_ptr;
typedef boost::variate_generator<boost::mt19937&, boost::uniform_real<> > unif_var;
typedef boost::shared_ptr<unif_var> unit_ptr;
typedef boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > norm_var;
typedef boost::shared_ptr<norm_var> norm_ptr;
typedef boost::variate_generator<boost::mt19937&, boost::random::uniform_int_distribution<> > unid_var;


//
// Global variable for Particle
//
unsigned multistep = 0;
string runtag, outdir;
char threading_on = 0;
int myid = 0;
pthread_mutex_t mem_lock;

//
// Physical constants
//
double a0    = 5.2917721092e-9;	    // cm (2xBohr radius)
double boltz = 1.3806488e-16;	    // cgs
double mp    = 1.67262178e-24;	    // g
double amu   = 1.660011e-24;	    // atomic mass unit
double pc    = 3.08567758e18;	    // Parsec (in cm)
double msun  = 1.9891e33;	    // Solar mass (in g)
double year  = 365.242*24.0*3600.0; // seconds per year

std::string atomic_specie[3] = {"electron",     "H",     "He"    };
double      atomic_masses[3] = {0.000548579909, 1.00794, 4.002602};

double Lunit;
double Tunit;
double Vunit;
double Munit;

std::vector<double> LL;

//
// boost RNGs
//
gen_ptr     gen;
uniform_ptr uniform;
normal_ptr  normal;
unit_ptr    Unit;
norm_ptr    Norm;

// ION collide types
//
enum Itype { Hybrid, Trace, Weight, Direct };
std::map<std::string, Itype> Types
{ {"Hybrid", Hybrid}, {"Trace", Trace}, {"Weight", Weight}, {"Direct", Direct} };

// Use CHIANTI or ION for ionization-recombination equilibrium
//
bool use_chianti = false;

/**
   Make Uniform temperature box of gas
*/
void InitializeUniform(std::vector<Particle>& p, double mass, double molW,
                       std::map<unsigned char, double>& T, vector<double> &L,
		       Itype type, int ne)
{
  unsigned npart = p.size();
  double   rho   = mass/(L[0]*L[1]*L[2]);
  
  std::cout << std::string(70, '-')                 << std::endl;
  if (T.size()>1 and type != Trace) {
    for (auto v : T) {
      std::ostringstream sout;
      sout << "Temp " << atomic_specie[v.first] << ":";
      std::cout << std::left << std::setw(13) << sout.str() << v.second
		<< " K "   << std::endl;
    }
  } else {
    std::cout << "Temp:        "  << T[0]  << " K "   << std::endl;
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
  std::map<unsigned char, double> varI, varE;
  for (auto v : T) {
    unsigned char Z = v.first;
    if (Z>0) {
      varI[Z] = sqrt((boltz*T[Z])/(atomic_masses[Z]*amu)) / Vunit; // Ion
      varE[Z] = sqrt((boltz*T[Z])/(atomic_masses[0]*amu)) / Vunit; // Electron
    } else {
      varI[Z] = sqrt((boltz*T[0])/(molW*amu))             / Vunit; // Fiducial particle
      varE[Z] = sqrt((boltz*T[0])/(atomic_masses[0]*amu)) / Vunit; // Electrons
    }
  }
  
  double tKEi  = 0.0;
  double tKEe  = 0.0;
  double numbI = 0.0;
  double numbE = 0.0;
  double Eunit = Munit*Vunit*Vunit;
  double ttemp = T[0];
  
  size_t nd = 0;		// Get the species fraction location
				// in dattrib
  if (type == Hybrid) {
    std::ifstream hin ("species.spec");
    std::string method, line;
    std::getline (hin, method);
    std::getline (hin, line);
    std::istringstream sin(line);
    sin >> nd;
    sin >> nd;
  }

  for (unsigned i=0; i<npart; i++) {
    
    double KE = 0.0;
    
    if (type == Trace) {
      
      for (unsigned k=0; k<3; k++) {
	p[i].pos[k] = L[k]*(*Unit)();
	p[i].vel[k] = varI[0]*(*Norm)();
      }
      
      if (ne>=0) {
	for (int l=0; l<3; l++) {
	  p[i].dattrib[ne+l] = varE[0] * (*Norm)();
	}
      }

    } else {
      
      KeyConvert kc(p[i].iattrib[0]);
      speciesKey sKey  = kc.getKey();
      unsigned short Z = sKey.first;
      unsigned short C = sKey.second;
      
      for (unsigned k=0; k<3; k++) {
	p[i].pos[k] = L[k]*(*Unit)();
	p[i].vel[k] = varI[Z] * (*Norm)();
	KE += p[i].vel[k] * p[i].vel[k];
      }
      
      double KEe = 0.0;
      if (ne>=0) {
	for (int l=0; l<3; l++) {
	  p[i].dattrib[ne+l] = varE[Z] * (*Norm)();
	  KEe += p[i].dattrib[ne+l] * p[i].dattrib[ne+l];
	}
      }
      
      // Ion KE
      //
      tKEi  += 0.5 * p[i].mass * KE * Eunit;

      // Ion number
      //
      numbI += p[i].mass/atomic_masses[Z] * Munit / amu;

      // Electron KE and number
      //
      if (type == Hybrid) {
	double eta = 0.0;
	for (unsigned short C=1; C<=Z; C++)
	  eta += p[i].dattrib[nd+C]*C;
	eta *= p[i].mass/atomic_masses[Z];
	numbE += eta * Munit / amu;
	tKEe  += 0.5 * eta * atomic_masses[0] * KEe * Eunit;
      } else {
	tKEe  += 0.5 * p[i].mass * atomic_masses[0]/atomic_masses[Z] * KEe * Eunit;
	numbE += p[i].mass/atomic_masses[Z] * Munit / amu * C;
      }

      ttemp = T[Z];
    }
    
    if (p[i].dattrib.size()>0) p[i].dattrib[0] = ttemp;
    else p[i].dattrib.push_back(ttemp);
    
    if (p[i].dattrib.size()>1) p[i].dattrib[1] = rho;
    else p[i].dattrib.push_back(rho);
    
    if (type == Direct) {
      if (p[i].dattrib.size()>2) p[i].dattrib[2] = KE;
      else p[i].dattrib.push_back(KE);
    }
  }
  
  if (type == Hybrid) {
    std::cout << "T (ion):      " << tKEi/(1.5*numbI*boltz) << std::endl
	      << "T (elec):     " << tKEe/(1.5*numbE*boltz) << std::endl
	      << "N (elec/ion): " << numbE/numbI            << std::endl
	      << std::string(70, '-') << std::endl;
  }
  else if (type != Trace) {
    std::cout << "T (ion):     " << tKEi/(1.5*numbI*boltz) << std::endl
	      << "T (elec):    " << tKEe/(1.5*numbI*boltz) << std::endl
	      << std::string(70, '-') << std::endl;
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
  
  for (int n=0; n<particles.size(); n++) {
    
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
		  << std::setw(16) << std::get<0>(i.second)*Munit/(atomic_masses[i.first.first]*mp*vol)
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
 double M, std::map<unsigned char, double>& T, int ne, int ni, int nd)
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
	   << " -T " << T[n] << " -o " << ioneq;
      
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
	   << " -T " << T[n];
      
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
  double norm = std::accumulate(sF.begin(), sF.end(), 0.0);
  if (fabs(norm - 1.0)>1.0e-16) {
    std::cout << "Normalization change: " << norm << std::endl;
  }
  
  std::vector<double> frcS(NS), cumS(NS);
  for (size_t i=0; i<NS; i++) {
    sF[i]  /= norm;
    frcS[i] = sF[i]/atomic_masses[sZ[i]];
    cumS[i] = frcS[i] + (i ? cumS[i-1] : 0);
  }
  
  double normC = cumS.back();
  
  for (size_t i=0; i<NS; i++) {
    frcS[i] /= normC;
    cumS[i] /= normC;
  }
  
  for (size_t i=0; i<N; i++) {
    
    double rz = (*Unit)();
    double rc = (*Unit)();
    
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
    for (size_t j=0; j<Zi+1; j++) {
      if (rc < cuml[indx][j]) {
	Ci = j+1;
	break;
      }
    }
    
    particles[i].mass  = M/N * atomic_masses[sZ[indx]] * normC;
    
    particles[i].iattrib.resize(ni, 0);
    particles[i].dattrib.resize(nd, 0);
    
    if (ne>=0) {			 // Add the use_elec fields
      for (int l=0; l<4; l++) particles[i].dattrib.push_back(0.0);
    }
    
    KeyConvert kc(speciesKey(Zi, Ci));
    particles[i].iattrib[0] = kc.getInt();
  }
  
  std::ofstream out("species.spec");
  out << "direct" << std::endl;
  out << std::setw(6) << nd << std::endl;
  for (size_t indx=0; indx<NS; indx++) { 
    out << std::setw(6) << static_cast<unsigned>(sZ[indx])
	<< std::endl;
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
 double M, std::map<unsigned char, double>& T, int& ne, int ni, int nd)
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
	   << " -T " << T[n] << " -o n" << ioneq;
      
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
	   << " -T " << T[n];
      
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
    wght[i] = frcS[i]/atomic_masses[sZ[i]];
  
  boost::random::uniform_int_distribution<> dist(0, NS-1);
  unid_var unifd(*gen, dist);
  
  for (size_t i=0; i<N; i++) {
    // Get the species
    size_t indx = unifd();
    unsigned char Ci = 1, Zi = sZ[indx];
    double rc = (*Unit)();
    
    // Get the ionization state
    for (size_t j=0; j<Zi+1; j++) {
      if (rc < cuml[indx][j]) {
	Ci = j+1;
	break;
      }
    }
    
    particles[i].mass  = M/N * sF[indx] * NS;
    
    particles[i].iattrib.resize(ni, 0);
    particles[i].dattrib.resize(nd, 0);
    particles[i].dattrib.push_back(0.0); // Add the use_cons field
    if (ne>=0) {			 // Add the use_elec fields
      for (int l=0; l<4; l++) particles[i].dattrib.push_back(0.0);
    }
    
    KeyConvert kc(speciesKey(Zi, Ci));
    particles[i].iattrib[0] = kc.getInt();
  }
  
  std::ofstream out("species.spec");
  
  out << "weight" << std::endl;
  out << std::setw(6) << nd++;
  if (ne>=0) out << std::setw(6) << (ne=nd);
  out << std::endl;
  
  for (size_t indx=0; indx<NS; indx++) { 
    out << std::setw(6)  << static_cast<unsigned>(sZ[indx])
	<< std::setw(16) << wght[indx]
	<< std::setw(16) << M/N * sF[indx] * NS
	<< std::endl;
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
 std::vector< std::vector<double> >& sI,
 double M, std::map<unsigned char, double>& T, int& ne, int ni, int nd)
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
	   << " -T " << T[n] << " -o " << ioneq;
      
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
	  cuml.push_back(v);
	}
      }
    }
    else {
      
      std::ostringstream sout;
      sout << "mpirun -np 1 genIonRecomb"
	   << " -Z " << static_cast<unsigned>(n)
	   << " -T " << T[n];
      
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
  
  auto it = std::max_element(std::begin(sZ), std::end(sZ));
  size_t maxSp = *it;
  
  wght[0] = W_H;
  for (size_t i=1; i<NS; i++)
    wght[i] = frcS[i]/atomic_masses[sZ[i]];
  
  boost::random::uniform_int_distribution<> dist(0, NS-1);
  unid_var unifd(*gen, dist);
  
  for (size_t i=0; i<N; i++) {
    // Get the species
    size_t indx = unifd();
    unsigned char Ci = 0, Zi = sZ[indx];
    
    particles[i].mass  = M/N * sF[indx] * NS;
    
    particles[i].iattrib.resize(ni, 0);
    particles[i].dattrib.resize(nd, 0);

    // Add the use_cons field
    particles[i].dattrib.push_back(0.0);

    // Add the ionization states
    for (auto v : frac[Zi]) {
      particles[i].dattrib.push_back(v);
    }

    // Pad
    for (size_t v=frac[Zi].size(); v<=maxSp; v++) 
      particles[i].dattrib.push_back(0.0);

    // Add the use_elec fields
    if (ne>=0) {
      for (int l=0; l<4; l++) particles[i].dattrib.push_back(0.0);
    }
    
    KeyConvert kc(speciesKey(Zi, Ci));
    particles[i].iattrib[0] = kc.getInt();
  }
  
  std::ofstream out("species.spec");
  
  out << "hybrid" << std::endl;
  out << std::setw(6) << nd++;
  out << std::setw(6) << nd; nd += maxSp + 1;
  if (ne>=0) out << std::setw(6) << (ne=nd);
  out << std::endl;
  
  for (size_t indx=0; indx<NS; indx++) { 
    out << std::setw(6)  << static_cast<unsigned>(sZ[indx])
	<< std::setw(16) << wght[indx]
	<< std::setw(16) << M/N * sF[indx] * NS
	<< std::endl;
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
 std::vector<double>& sF, double M, double T,
 int& ne, int ni, int nd)
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
	   << " -T " << T << " -o " << ioneq;
      
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
	   << " -T " << T;
      
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
      }
    }
  }
  
  int N = particles.size();
  
  // Compute cumulative species
  // distribution
  size_t NS = sF.size();
  // Normalize sF
  double norm = std::accumulate(sF.begin(), sF.end(), 0.0);
  if (fabs(norm - 1.0)>1.0e-16) {
    std::cout << "Normalization change: " << norm << std::endl;
  }
  for (size_t i=0; i<NS; i++) sF[i]  /= norm;
  
  for (size_t indx=0; indx<NS; indx++) { 
    double norm = std::accumulate(frac[indx].begin(), frac[indx].end(), 0.0);
    for (auto &i : frac[indx]) i /= norm;
  }
  
  for (size_t i=0; i<N; i++) {
    
    particles[i].mass  = M/N;
    
    particles[i].iattrib.resize(ni, 0);
    particles[i].dattrib.resize(nd, 0);
    
    // Get the species
    double test = 0.0;
    for (size_t indx=0; indx<NS; indx++) { 
      // Get the ionization state
      for (size_t j=0; j<sZ[indx]+1; j++) {
	particles[i].dattrib.push_back(sF[indx]*frac[indx][j]);
	test += sF[indx] * frac[indx][j];
      }
    }
    assert ( fabs(test-1.0) < 1.0e-12 );

    // Add the use_elec fields
    if (ne>=0) {
      for (int l=0; l<4; l++) particles[i].dattrib.push_back(0.0);
    }
  }
  
  if (ne>=0) {
    ne = nd;
    for (size_t indx=0; indx<NS; indx++) { 
      for (size_t j=0; j<sZ[indx]+1; j++) ne++;
    }
  }

  std::ofstream out("species.spec");
  out << "trace" << std::endl;
  // Electron position (-1 for none)
  out << std::setw(6) << ne << std::endl;
  int cntr = nd;
  for (size_t indx=0; indx<NS; indx++) { 
    for (size_t j=0; j<sZ[indx]+1; j++) {
      out << std::setw(6) << static_cast<unsigned>(sZ[indx])
	  << std::setw(6) << j + 1
	  << std::setw(6) << cntr++
	  << std::endl;
    }
  }
}

int main (int ac, char **av)
{
  Itype type = Direct;
  double D, L, temp;
  std::string config;
  std::string oname;
  unsigned seed;
  int ne = -1, ni = 2, nd = 6;
  int npart;
  
  std::string cmd_line;
  for (int i=0; i<ac; i++) {
    cmd_line += av[i];
    cmd_line += " ";
  }
  
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("electrons",       "set up for weighted or hybrid species with electrons")
    ("CHIANTI,C",	po::value<bool>(&use_chianti)->default_value(false),
     "use CHIANTI to set recombination-ionization equilibrium")
    ("temp,T",		po::value<double>(&temp)->default_value(25000.0),
     "set the temperature in K")
    ("dens,D",		po::value<double>(&D)->default_value(1.0),
     "density in particles per cc")
    ("length,L",	po::value<double>(&L)->default_value(1.0),
     "length in system units")
    ("number,N",	po::value<int>(&npart)->default_value(250000),
     "number of particles")
    ("seed,s",		po::value<unsigned>(&seed)->default_value(11),
     "random number seed")
    ("Lunit,l",		po::value<double>(&Lunit)->default_value(1.0),
     "length scale of the system in pc")
    ("Tunit,t",		po::value<double>(&Tunit)->default_value(1.0e5),
     "time scale in years")
    ("Munit,m",		po::value<double>(&Munit)->default_value(0.1),
     "mass scale in solar masses")
    ("num-int,i",	po::value<int>(&ni)->default_value(2),
     "number of integer attributes")
    ("num-double,d",	po::value<int>(&nd)->default_value(6),
     "number of double attributes")
    ("config,c",	po::value<std::string>(&config)->default_value("makeIon.config"),
     "element config file")
    ("output,o",	po::value<std::string>(&oname)->default_value("out.bods"),
     "body output file")
    ;
  
  
  po::variables_map vm;
  
  try {
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);    
  } catch (po::error& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }
  
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    std::cout << "Example: " << std::endl;
    std::cout << "\t" << av[0]
	      << " --temp=25000 --number=250000 --output=out.bod" << std::endl;
    return 1;
  }
  
  if (vm.count("hybrid")) {
    type = Hybrid;
  }
  
  if (vm.count("weight")) {
    type = Weight;
  }
  
  if (vm.count("trace")) {
    type = Trace;
  }
  
  if (vm.count("electrons")) {
    ne   = 0;
  }
  
  if (myid==0) {
    std::string prefix("makeIon");
    std::string cmdFile = prefix + ".cmd_line";
    std::ofstream out(cmdFile.c_str());
    if (!out) {
      std::cerr << "makeIon: error opening <" << cmdFile
                << "> for writing" << std::endl;
    } else {
      out << cmd_line << std::endl;
    }
  }
  
  Lunit  *= pc;
  Tunit  *= year;
  Munit  *= msun;
  Vunit   = Lunit/Tunit;
  
  gen     = gen_ptr    (new boost::mt19937(seed));
  uniform = uniform_ptr(new boost::uniform_real<>(0.0, 1.0));
  normal  = normal_ptr (new boost::normal_distribution<>(0.0, 1.0));
  Unit    = unit_ptr   (new unif_var(*gen, *uniform));
  Norm    = norm_ptr   (new norm_var(*gen, *normal));
  
  //
  // Define the atomic species statistics
  //
  
  // Short alias for this namespace
  namespace pt = boost::property_tree;
  
  // Element data
  std::vector<unsigned char> sZ;
  std::vector<double>        sF;

  // Temperatures
  std::map<unsigned char, double> T;

  // Parse element file
  {
    // Create a root
    pt::ptree iroot;
  
    // Load the json file in this ptree
    pt::read_json(config, iroot);
  
    double norm = 0.0;

    std::string st = iroot.get("type", "Direct");
    if (Types.find(st) == Types.end()) {
      std::cout << "Type <" << st << "> is not valid" << std::endl
		<< "Valid types are:";
      for (auto v : Types) std::cout << " " << v.first;
      std::cout << std::endl;
      exit(-1);
    }

    type = Types[st];

    if (type==Trace) {
      T[0] = iroot.get("temp", 100000.0);
    }

    for (pt::ptree::value_type &row : iroot.get_child("elements")) {
    
      std::istringstream sin(row.first);
      unsigned i;
      sin >> i;
    
      sZ.push_back(i);
      sF.push_back(row.second.get<double>("mfrac"));
      norm += sF.back();
      if (type!=Trace) T[i] = row.second.get<double>("temp");
    }

    if (norm>0.0) {
      for (auto & v : sF) v /= norm;
    } else {
      std::cout << "Error: zero mass fraction norm" << std::endl;
      exit(-1);
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
  double vol = 1.0;
  for (auto v : LL) vol *= v*Lunit;
  
  // Mass in box in m_p
  //
  double Mass = mp*D*vol/Munit;
  
  vector<Particle> particles(npart);
  
  double molW = 1.0;

  // Initialize the Z, C's	
  //
  switch (type) {
  case Hybrid:
    InitializeSpeciesHybrid(particles, sZ, sF, sI, Mass, T, ne, ni, nd);
    break;
  case Weight:
    InitializeSpeciesWeight(particles, sZ, sF, sI, Mass, T, ne, ni, nd);
    break;
  case Trace:
    InitializeSpeciesTrace (particles, sZ, sF, Mass, temp, ne, ni, nd);
    // Compute molecular weight
    molW = 0.0;
    for (size_t k=0; k<sZ.size(); k++) molW += sF[k]/atomic_masses[sZ[k]];
    molW = 1.0/molW;
    break;
  case Direct:
  default:
    InitializeSpeciesDirect(particles, sZ, sF, sI, Mass, T, ne, ni, nd);
  }
  
  // Initialize the phase space vector
  //
  InitializeUniform(particles, Mass, molW, T, LL, type, ne);
  
  // Output file
  //
  writeParticles(particles, oname, type, sF, sI);
  
  return 0;
}
