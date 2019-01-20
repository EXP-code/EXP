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

#include <yaml-cpp/yaml.h>

#include "Particle.H"
#include "globalInit.H"
#include "Species.H"
#include "atomic_constants.H"

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

double Lunit;
double Tunit;
double Vunit;
double Munit;

std::vector<double> LL;

PeriodicTable PT;

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
bool use_chianti   = false;
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
void InitializeUniform(std::vector<Particle>& p, std::vector<double>& mass, double molW,
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

    For mean-mass trace-species algorithm:

    3/2*N_e*k_B*T = 3/2*m*k_B*T/mu(Z) = 3/2*m/mu(Z)*eta*m_e*v_e^2

    or 
    
    v_e^2 = k_B*T/(m_e*eta)
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
      // unsigned short C = sKey.second;
      
      for (unsigned k=0; k<3; k++) {
	p[i].pos[k] = L[k]*(*Unit)();
	p[i].vel[k] = varI[Z] * (*Norm)();
	KE += p[i].vel[k] * p[i].vel[k];
      }
      
      if (ne>=0) {
	for (int l=0; l<3; l++) {
	  p[i].dattrib[ne+l] = varE[Z] * (*Norm)();
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

  boost::random::uniform_int_distribution<> dist(0, 1);
  unid_var selP(*gen, dist);

  for (unsigned i=0; i<npart; i++) {
    
    double KE = 0.0;
    
    unsigned char wh = static_cast<unsigned char>(selP());


    p[i].iattrib.resize(ni, 0);
    p[i].dattrib.resize(nd, 0);

    if (type == Trace) {
      
      for (unsigned k=0; k<3; k++) {
	if (k==0)
	  p[i].pos[k] = 0.5*L[k]*((*Unit)() + wh);
	else
	  p[i].pos[k] = L[k]*(*Unit)();
	p[i].vel[k] = varI[wh][0]*(*Norm)();
      }
      
      if (ne>=0) {
	for (int l=0; l<3; l++) {
	  p[i].dattrib[ne+l] = varE[wh][0] * (*Norm)();
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
	  p[i].pos[k] = 0.5*L[k]*((*Unit)() + wh);
	else
	  p[i].pos[k] = L[k]*(*Unit)();
	p[i].vel[k] = varI[wh][Z] * (*Norm)();
	KE += p[i].vel[k] * p[i].vel[k];
      }
      
      if (ne>=0) {
	for (int l=0; l<3; l++) {
	  p[i].dattrib[ne+l] = varE[wh][Z] * (*Norm)();
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
  
  boost::random::uniform_int_distribution<> dist(0, NS-1);
  unid_var unifd(*gen, dist);
  
  for (int i=0; i<N; i++) {
    // Get the species
    size_t indx = unifd();
    unsigned char Ci = 1, Zi = sZ[indx];
    double rc = (*Unit)();
    
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
  
  boost::random::uniform_int_distribution<> dist(0, NS-1);
  unid_var unifd(*gen, dist);
  
  for (int i=0; i<N; i++) {
    // Get the species
    size_t indx = unifd();
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
 std::vector<double>& sF, std::vector<double>& M,
 std::vector< std::map<unsigned char, double> >& T,
 std::vector<double>& L,
 Mtype model, int sp, int ne, bool mm)
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
	for (int C=0; C<n; C++) {
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
  
  // Compute cumulative species
  // distribution
  int NS = sF.size();

  double molW = 0.0;
  for (int i=0; i<NS; i++)
    molW += sF[i]/PT[sZ[i]]->weight();
  molW = 1.0/molW;
  
  // Normalize sF
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

  // Setup for mean-mass correction
  //
  if (mm) {
    std::cout << std::string(70, '-') << std::endl;

    for (size_t nc=0; nc<Ncomp; nc++) {
      eta[nc] = 0.0;
      double nrm = 0.0;
      for (int indx=0; indx<NS; indx++) { 
	int C = 0;
	for (auto v : frac[nc][indx])  {
	  double wgt = sF[indx]/PT[sZ[indx]]->weight() * v;
	  eta[nc] += wgt * C++;
	  nrm     += wgt;
	}
      }
      if (nrm>0.0) eta[nc] /= nrm;
      std::ostringstream lab; lab << "Eta (" << nc << "):";
      std::cout << std::left << std::setw(13) << lab.str()
		<< eta[0] << std::endl;
      lab.str(""); lab << "Mol (" << nc << "):";
      std::cout << std::left << std::setw(13) << lab.str()
		<< 1.0/nrm << std::endl;
    }
  }
  
  double tKEi = 0.0, tKEe = 0.0, numb = 0.0;

  for (int i=0; i<N; i++) {
    
    size_t wh = 0;
    if (model == Interface) {
      if (particles[i].pos[0] > 0.5*L[0]) wh = 1;
    }

    particles[i].mass  = M[wh]/N;
    
    // Sanity check
    double test = 0.0;

    int cur = sp;
    // Get the species
    //
    for (int indx=0; indx<NS; indx++) { 
      // Get the ionization state
      //
      for (int j=0; j<sZ[indx]+1; j++) {
	particles[i].dattrib[cur++] = sF[indx]*frac[wh][indx][j];
	test += sF[indx] * frac[wh][indx][j];
      }
    }
    assert ( fabs(test-1.0) < 1.0e-12 );

    double KEi = 0.0, KEe = 0.0;
    for (int k=0; k<3; k++) {
      KEi += particles[i].vel[k] * particles[i].vel[k];
      if (ne>=0) {
	particles[i].dattrib[ne+k] /= sqrt(eta[wh]); // mean-mass correction
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
      for (int j=0; j<sZ[indx]+1; j++) {
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
      for (int j=0; j<sZ[indx]+1; j++) {
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
  Itype type  = Direct;
  Mtype model = Uniform;
  double D, L, R, Temp, Telec;
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
  
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("electrons",       "set up for weighted or hybrid species with electrons")
    ("meanmass",        "set up for the mean-mass algorithm")
    ("yaml",            "write YAML species config file")
    ("old",             "write old-style species config file")
    ("CHIANTI,C",	po::value<bool>(&use_chianti)->default_value(false),
     "use CHIANTI to set recombination-ionization equilibrium")
    ("INIT,I",	        po::value<bool>(&use_init_file)->default_value(false),
     "use init file to set recombination-ionization equilibrium")
    ("dens,D",		po::value<double>(&D)->default_value(1.0),
     "density in particles per cc")
    ("temp,T",		po::value<double>(&Temp)->default_value(-1.0),
     "override config file temperature for Trace, if >0")
    ("Telec",		po::value<double>(&Telec)->default_value(-1.0),
     "temperature for electrons, if Telec>0")
    ("length,L",	po::value<double>(&L)->default_value(1.0),
     "length in system units")
    ("ratio,R",		po::value<double>(&R)->default_value(1.0),
     "slab length ratio (1 is cube")
    ("number,N",	po::value<int>(&npart)->default_value(250000),
     "number of particles")
    ("seed,s",		po::value<unsigned>(&seed)->default_value(11),
     "random number seed")
    ("Lunit,l",		po::value<double>(&Lunit)->default_value(1.0),
     "length scale of the system in pc")
    ("Tunit,t",		po::value<double>(&Tunit)->default_value(1.0e3),
     "time scale in years")
    ("Munit,m",		po::value<double>(&Munit)->default_value(1.0),
     "mass scale in solar masses")
    ("num-int,i",	po::value<int>(&ni)->default_value(2),
     "number of integer attributes")
    ("num-double,d",	po::value<int>(&nd)->default_value(6),
     "base number of double attributes")
    ("config,c",	po::value<std::string>(&config)->default_value("makeIon.config"),
     "element config file")
    ("output,o",	po::value<std::string>(&oname)->default_value("out"),
     "output prefix")
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
	      << "--length=0.0001 --number=1000 --dens=0.01 --electrons -c trace.config --output=out"
	      << std::endl;
    return 1;
  }
  
  bool mm = false;
  if (vm.count("meanmass")) {
    mm = true;
  }

  if (vm.count("yaml")) {
    use_yaml = true;
  }

  if (vm.count("old")) {
    use_yaml = false;
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
  std::vector< std::map<unsigned char, double> > T;
  std::vector<double> rho;

  // Species position
  int sp = -1;

  // Default molecular weight
  double molW = 1.0;

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

    std::string md = iroot.get("model", "Uniform");
    if (Models.find(md) == Models.end()) {
      std::cout << "Model <" << md << "> is not valid" << std::endl
		<< "Valid types are:";
      for (auto v : Models) std::cout << " " << v.first;
      std::cout << std::endl;
      exit(-1);
    }

    model = Models[md];

    if (iroot.get("electrons", true)) {
      ne = 0;
    }
  
    if (model==Interface) {
      T.resize(2);
    } else {
      T.resize(1);
    }

    double fH = 0.0;

    for (pt::ptree::value_type &row : iroot.get_child("elements")) {
    
      std::istringstream sin(row.first);
      unsigned i;
      sin >> i;
    
      sZ.push_back(i);

      if (row.second.count("logN")) {
	if (i==1) fH = row.second.get<double>("logN");
	double wgt = pow(10.0, row.second.get<double>("logN") - fH)*PT[i]->weight();
	sF.push_back(wgt);
	norm += sF.back();
      } else if (row.second.count("mfrac")) {
	sF.push_back(row.second.get<double>("mfrac"));
	norm += sF.back();
      } else {
	std::cerr << "Missing element definition for atomic number " << i
		  << std::endl;
	exit(-3);
      }
      if (type!=Trace) T[0][i] = row.second.get<double>("temp");
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
      try {
	rho[0]  = iroot.get<double>("components.1.Density");
	rho[1]  = iroot.get<double>("components.2.Density");
	T[0][0] = iroot.get<double>("components.1.Temp");
	T[0][1] = iroot.get<double>("components.1.Etemp", T[0][0]);
	T[1][0] = iroot.get<double>("components.2.Temp");
	T[1][1] = iroot.get<double>("components.2.Etemp", T[1][0]);
      }
      catch (pt::ptree_error & error) {
	std::cerr << "Error: " << error.what() << std::endl;
	exit(-2);
      }
    } else {
      rho.push_back(D);
      if (Temp>0.0) T[0][0]  = Temp;
      else          T[0][0]  = iroot.get("temp", 100000.0);
      T[0][1] = T[0][0];
      if (Telec>0.0) T[0][1] = Telec;
      else           T[0][1] = iroot.get("telc", T[0][0]);
    }
      
    if (type==Trace) {
      nd++;			// Conservation of energy position
      sp = nd;

				// Augment for species number
      for (size_t indx=0; indx<sF.size(); indx++) { 
	for (int j=0; j<sZ[indx]+1; j++) nd++;
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
  
  vector<Particle> particles(npart);
  
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
    InitializeUniform(particles, Mass, molW, T, LL, type, sp, ne, ni, nd);
    break;
  }

  // Initialize the Z, C's	
  //
  switch (type) {
  case Hybrid:
    InitializeSpeciesHybrid(particles, sZ, sF, sI, Mass, T, sp, ne);
    break;
  case Weight:
    InitializeSpeciesWeight(particles, sZ, sF, sI, Mass, T, sp, ne);
    break;
  case Trace:
    InitializeSpeciesTrace (particles, sZ, sF, Mass, T, LL, model, sp, ne, mm);
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
