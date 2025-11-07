/*
  A tapered Mestel disk IC generator
*/
                                // C++/STL headers
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include <string>
#include <vector>
#include <array>

#include <omp.h>

#include "mestel.H"

#include "Progress.H"		// Progress bar
#include "cxxopts.H"		// Option parsing

int 
main(int ac, char **av)
{
  //=====================
  // Begin option parsing
  //=====================

  int          N;		// Number of particles
  int          Nrepl;		// Number of particle replicates per orbit
  double       mu, nu, Ri, Ro;	// Taper paramters
  double       Rmin, Rmax;      // Radial range
  double       sigma;		// Velocity dispersion
  std::string  bodyfile;	// Output file
  unsigned     seed;		// Will be inialized by /dev/random if
				// not set on the command line

  cxxopts::Options options(av[0], "Ideal tapered Mestel IC generator");

  options.add_options()
    ("h,help",    "Print this help message")
    ("V,nozerovel", "Do not zero the mean velocity")
    ("P,nozeropos", "Do not zero the center of mass")
    ("d,debug",   "Print debug grid")
    ("N,number",  "Number of particles to generate",
     cxxopts::value<int>(N)->default_value("100000"))
    ("n,nu",      "Inner taper exponent (0 for no taper)",
     cxxopts::value<double>(nu)->default_value("2.0"))
    ("m,mu",      "Outer taper exponent (0 for no taper)",
     cxxopts::value<double>(mu)->default_value("2.0"))
    ("i,Ri",      "Inner radius for taper",
     cxxopts::value<double>(Ri)->default_value("1.0"))
    ("o,Ro",      "Outer radius for taper",
     cxxopts::value<double>(Ro)->default_value("20.0"))
    ("r,Rmin",    "Inner radius for model",
     cxxopts::value<double>(Rmin)->default_value("0.001"))
    ("R,Rmax",    "Outer radius for model",
     cxxopts::value<double>(Rmax)->default_value("50.0"))
    ("S,sigma",   "Radial velocity dispersion",
     cxxopts::value<double>(sigma)->default_value("1.0"))
    ("s,seed",    "Random number seed. Default: use /dev/random",
     cxxopts::value<unsigned>(seed))
    ("q,Nrepl",   "Number of particle replicates per orbit",
     cxxopts::value<int>(Nrepl)->default_value("1"))
    ("f,file",    "Output body file",
     cxxopts::value<std::string>(bodyfile)->default_value("zang.bods"))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(ac, av);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    std::cout << options.help() << std::endl << std::endl;
    return 1;
  }

  // Set from /dev/random if not specified
  if (vm.count("seed")==0) {
    seed = std::random_device{}();
  }

  // Set particle number consitent with even replicants
  if (Nrepl<1) Nrepl = 1;
  if (Nrepl>1) N = (N/Nrepl)*Nrepl;

  // Make sure N>0
  if (N<=0) {
    std::cerr << av[0] << ": you must requiest at least one body"
	      << std::endl;
  }

  // Open the output file
  //
  std::ofstream out(bodyfile);
  if (not out) {
    std::string msg(av[0]);
    msg +=  ": output file <" + bodyfile + "> can not be opened";
    throw std::runtime_error(msg);
  }

  SphericalOrbit::ZFRAC=0.3;	// TEST

  // Create the model
  //
  auto model = std::make_shared<TaperedMestelDisk>(nu, mu, Ri, Ro,
						   1.0, Rmin, Rmax);
  model->setup_df(sigma);

  // Replicant logic
  //
  int Number = N/Nrepl;
  double dPhi = 2.0*M_PI/Nrepl;

  // Progress bar
  //
  std::shared_ptr<progress::progress_display> progress;
  int nomp = 1;
#pragma omp parallel
  {
    nomp = omp_get_num_threads();
    if (omp_get_thread_num()==0) {
      progress = std::make_shared<progress::progress_display>(Number/nomp);
    }
  }

  // Create an orbit grid
  //
  std::vector<std::shared_ptr<SphericalOrbit>> orb(nomp);
  for (auto & v : orb) v = std::make_shared<SphericalOrbit>(model);

  double Ktol = 0.01;
  double Kmin = Ktol, Kmax = 1.0 - Ktol;

  double Emin = 0.5*Rmin*model->get_dpot(Rmin) + model->get_pot(Rmin);
  // double Emax = 0.5*Rmax*model->get_dpot(Rmax) + model->get_pot(Rmax);
  double Emax = model->get_pot(Rmax);

  // Scan to find the peak df
  //
  const int numE = 800;
  const int numK = 40;
  Eigen::VectorXd cumE(numE+1), cumF(numE+1), topF(numE+1);
  double peak = 0.0;
  double dE = (Emax - Emin)/numE, dK = (1.0 - 2.0*Ktol)/numK;
  for (int i=0; i<=numE; i++) {
    double E = Emin + dE*i;
    cumE(i) = E;
    cumF(i) = 0.0;
    topF(i) = 0.0;
    for (int j=0; j<=numK; j++) {
      double K = Kmin + dK*j;
      orb[0]->new_orbit(E, K);
      double F = model->distf(E, orb[0]->get_action(1)) / orb[0]->get_freq(0);
      peak = std::max<double>(peak, F);
      topF(i) = std::max<double>(topF(i), F);
      cumF(i) += F * orb[0]->Jmax()/orb[0]->get_freq(0);
    }
  }

  if (peak <= 0.0) {
    throw std::runtime_error(std::string(av[0]) + ": peak DF is zero!");
  }

  // Improve the acceptance rejection by computing the cumulative
  // distribution
  //
  for (int i=1; i<=numE; i++) cumF[i] += cumF[i-1];

  if (cumF[numE] <= 0.0) {
    throw std::runtime_error(std::string(av[0]) + ": no mass on cum DF grid!");
  }

  for (int i=0; i<=numE; i++) cumF[i] /= cumF[numE];

  Linear1d Ecum(cumF, cumE), Ftop(cumE, topF);

  std::mt19937 gen(seed);
  std::uniform_real_distribution<> uniform(0.0, 1.0);

  // Save the position and velocity vectors
  //
  std::vector<std::array<double, 3>> pos(N), vel(N);

  // Maximum number rejection-method iterations
  int itmax = 100000;

  // Track number of iteration overflows
  int over  = 0;

  // Position and velocity zeroing
  //
  std::vector<std::array<double, 3>> zeropos(nomp), zerovel(nomp);

  // Generation loop with OpenMP
  //
#pragma omp parallel for reduction(+:over)
  for (int n=0; n<Number; n++) {
    // Thread id
    int tid = omp_get_thread_num();

    // Loop variables
    double E, F, K;
    int j;
    for (j=0; j<itmax; j++) {

      E = Ecum.eval(uniform(gen));
      F = Ftop.eval(E);
      K = Kmin + (Kmax - Kmin)*uniform(gen);

      orb[tid]->new_orbit(E, K);
      double F = model->distf(E, orb[tid]->get_action(1)) / orb[tid]->get_freq(0);
      if (F/peak > uniform(gen)) break;
    }

    if (j==itmax) over++;

    double J   = orb[tid]->get_action(1);
    double T   = 2.0*M_PI/orb[tid]->get_freq(0)*uniform(gen);
    double r   = orb[tid]->get_angle(6, T);
    double w1  = orb[tid]->get_angle(1, T);
    double phi = 2.0*M_PI*uniform(gen) + orb[tid]->get_angle(7, T);

    double vt  = J/r;
    double vr  = sqrt(fabs(2.0*(E - model->get_pot(r)) - J*J/(r*r)));
    //                ^
    //                |
    // For sanity-----+

    if (w1 > M_PI) vr *= -1.0;	// Branch of radial motion

    for (int nn=0; nn<Nrepl; nn++) {
      int indx = n*Nrepl + nn;
      double Phi = phi + dPhi*nn;
      // Convert from polar to Cartesian
      //
      pos[indx][0] = r*cos(Phi);
      pos[indx][1] = r*sin(Phi);
      pos[indx][2] = 0.0;

      vel[indx][0] = vr*cos(Phi) - vt*sin(Phi);
      vel[indx][1] = vr*sin(Phi) + vt*cos(Phi);
      vel[indx][2] = 0.0;

      // Accumulate mean position and velocity
      //
      for (int k=0; k<3; k++) {
	zeropos[tid][k] += pos[indx][k];
	zerovel[tid][k] += vel[indx][k];
      }
    }
      
    // Print progress bar
    if (tid==0) ++(*progress);
  }
  std::cout << std::endl << "** Main loop complete" << std::endl;

  // Compute the particle mass
  //
  double mass = (model->get_mass(Rmax) - model->get_mass(Rmin))/N;

  // Reduce the mean position and velocity
  //
  for (int n=1; n<nomp; n++) {
    for (int k=0; k<3; k++) {
      zeropos[0][k] += zeropos[n][k];
      zerovel[0][k] += zerovel[n][k];
    }
  }
  for (int k=0; k<3; k++) {
    zeropos[0][k] /= N;
    zerovel[0][k] /= N;
  }

  std::cout << "** Position center: " << zeropos[0][0] << ", "
	    << zeropos[0][1] << ", " << zeropos[0][2] << std::endl;
  
  std::cout << "** Velocity center: " << zerovel[0][0] << ", "
	    << zerovel[0][1] << ", " << zerovel[0][2] << std::endl;
  
  if (vm.count("nozeropos"))
    for (int k=0; k<3; k++) zeropos[0][k] = 0.0;

  if (vm.count("nozerovel"))
    for (int k=0; k<3; k++) zerovel[0][k] = 0.0;

  std::cout << "** " << over << " particles failed acceptance" << std::endl
	    << "** Particle mass=" << mass << std::endl;

  out << std::setw(8) << N << std::setw(8) << 0 << std::setw(8) << 0
      << std::endl;

  double ektot = 0.0, clausius = 0.0;
  for (int n=0; n<N; n++) {
    out << std::setw(18) << mass;
    for (int k=0; k<3; k++) out << std::setw(18) << pos[n][k] - zeropos[0][k];
    for (int k=0; k<3; k++) out << std::setw(18) << vel[n][k] - zerovel[0][k];
    out << std::endl;

    double r2 = 0.0, v2 = 0.0;
    for (int k=0; k<3; k++) {
      r2 += pos[n][k]*pos[n][k];
      v2 += vel[n][k]*vel[n][k];
    }

    ektot += mass*v2;
    double r = sqrt(r2);
    clausius += mass*model->get_dpot(r)*r;
  }

  std::cout <<  "** 2T/VC=" << ektot/clausius << std::endl;

  if (vm.count("debug")) {
    std::cout << std::endl << "Peak per energy" << std::endl;
    for (int i=0; i<=numE; i++) {
      std::cout << std::setw(6)  << i
		<< std::setw(18) << cumE(i)
		<< std::setw(18) << cumF(i)
		<< std::setw(18) << topF(i)
		<< std::endl;
    }
    std::cout << std::endl;
  }


  return 0;
}
