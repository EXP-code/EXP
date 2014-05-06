/*
  Compute initial conditions for a uniform density gas
*/

#include <iostream>
#include <iomanip>

#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "Particle.H"
#include "Initialize.H"
#include "globalInit.H"

unsigned multistep = 0;

double ne          = 0.5;	      // electron density in n/cm^3
double pc          = 3.086e18;	      // Parsec (in cm)
double msun        = 1.9891e33;	      // Solar mass (in g)

double Lunit;
double Tunit;
double Vunit;
double Munit;

// boost RNGs

gen_ptr     gen;
uniform_ptr uniform;
normal_ptr  normal;
unit_ptr    Unit;
norm_ptr    Norm;

void writeParticles(std::vector<Particle>& particles, const string& file)
{
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

    out << setw(18) << particles[n].Z;
    out << setw(18) << particles[n].C;

    for (unsigned k=0; k<particles[n].iattrib.size(); k++)
      out << setw(12) << particles[n].iattrib[k];

    for (unsigned k=0; k<particles[n].dattrib.size(); k++)
      out << setw(18) << particles[n].dattrib[k];

    out << std::endl;
  }
}

void InitializeSpecies(const std::string& stateF,
		       std::vector<Particle> & particles, 
		       double T, 
		       std::vector<int> Zspecies, 
		       std::vector<double> fZ, 
		       int NZ)
{
  // double planck = 6.62606957e-27;
  double boltz = 1.381e-16;
  double me    = 9.10938e-28;
  double mp    = 1.67e-24;
  
  std::vector< std::vector< double> > cumRats;
  std::string inLine;
  std::cout << "Opening state file: " << stateF.c_str() << std::endl;
  std::ifstream sFile(stateF.c_str());
  std::cout << "Opened " << endl;
  if (sFile.is_open()) {
    while(sFile.good()) {
      std::vector<std::string> v;
      getline(sFile, inLine);
      std::cout << inLine << std::endl;
      std::istringstream iss(inLine);
      std::copy(std::istream_iterator<std::string>(iss), 
		std::istream_iterator<std::string>(), 
		std::back_inserter<std::vector<std::string> >(v));
      std::vector<double> ratZ;
      for (int i = 0; i < v.size(); i++) {
	ratZ.push_back(atof(v[i].c_str()));
	std::cout << v[i] << "\t";
      }
      std::cout << std::endl;
      cumRats.push_back(ratZ);
    }
  }

  // std::cout << cumRats.size() << endl;
  int N = particles.size();
  std::vector<double> fC(NZ);
  fC[0] = fZ[0];
  for(int i = 1; i < NZ; i++) {
    fC[i] = fC[i-1] + fZ[i];
  }
  for(int i = 0; i < N; i++) {
    double rz = (*Unit)();
    double rc = (*Unit)();
    int Ci = 1;
    int Zi;
    int indexZ;
    for (int k = 0; k < NZ; k++) { 
      if (rz < fC[k]) {
	Zi = Zspecies[k]; 
	indexZ = k; 
	k = NZ+5;
      }
    }
    for(int j = 0; j < Zi+1; j++) {
      if (rc < cumRats[indexZ][j]) {
	Ci = j+1;
	break;
      }
    }
    double mi = Zi*mp;
    particles[i].Z = Zi;
    particles[i].C = Ci;
  }
  
}

int main (int ac, char **av)
{
  double T, Lunit, Tunit, Munit;
  std::string stateF, oname;
  unsigned seed;
  int npart;

  po::options_description desc("Allowed options");
  desc.add_options()
    ("help,h",		"produce help message")
    ("temp,T",		po::value<double>(&T)->default_value(25000.0),
     "set the temperature in K")
    ("number,N",	po::value<int>(&npart)->default_value(250000),
     "set the temperature in K")
    ("seed,s",		po::value<unsigned>(&seed)->default_value(11),
     "random number seed")
    ("Lunit,l",		po::value<double>(&Lunit)->default_value(1.0),
     "length scale of the system in pc")
    ("Tunit,t",		po::value<double>(&Tunit)->default_value(1.0e5),
     "time scale in years")
    ("Munit,m",		po::value<double>(&Munit)->default_value(0.1),
     "mass scale in solar masses")
    ("ion,i",		po::value<std::string>(&stateF)->default_value("states.dat"),
     "ionization state file")
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
	      << " --temp=25000 --num=250000 --state=states.dat --out=uni_ISM.bod" << std::endl;
    return 1;
  }

  Lunit  *= pc;
  Tunit  *= 365.25*3600.0*24.0;
  Munit  *= msun;
  Vunit   = Lunit/Tunit;

  gen     = gen_ptr    (new boost::mt19937(seed));
  uniform = uniform_ptr(new boost::uniform_real<>(0.0, 1.0));
  normal  = normal_ptr (new boost::normal_distribution<>(0.0, 1.0));
  Unit    = unit_ptr   (new unif_var(*gen, *uniform));
  Norm    = norm_ptr   (new norm_var(*gen, *normal));


  //
  // Physical units
  //
  double a0    = 2.0*0.054e-7;	// cm (2xBohr radius)
  double boltz = 1.381e-16;	// cgs
  double mp    = 1.67e-24;	// g
  
  //
  // Define the atomic species statistics
  //
  const int Nspec = 2.0;
  std::vector< double > fspec(Nspec);
  std::vector< int > Zspec(Nspec);

  fspec[0] = 0.75; Zspec[0] = 1;
  fspec[1] = 0.25; Zspec[1] = 2;

  /* Other species
     fspec[2] = 0.0; Zspec[2] = 3; //Li
     fspec[3] = 0.0; Zspec[3] = 6; //C
     fspec[4] = 0.0; Zspec[4] = 7; //N
     fspec[5] = 0.0  Zspec[5] = 8; //O
     fspec[6] = 0.0; Zspec[6] = 12; //Mg
  */
  double D = 1.0;
  double P = 1.0;
  double L = 1.0;

  std::vector< double > LL(3);
  LL[0] = L; LL[1] = L; LL[2] = L;

  double Mass = D*(LL[0]*LL[1]*LL[2]);

  vector<Particle> particles(npart);

  // Initialize the Z, C's	

  InitializeSpecies(stateF, particles, T, Zspec, fspec, Nspec);
  
  InitializeUniform(particles, Mass, T, LL);
  
  writeParticles(particles, oname);

  return 0;
}
