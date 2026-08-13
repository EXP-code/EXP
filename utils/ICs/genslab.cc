// Generate slab initial conditions

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <random>
#include <string>
#include <memory>
#include <cmath>

#include <Eigen/Dense>

#include "massmodel1d.H"
#include "interp.H"
#include "cxxopts.H"

int
main(int argc, char **argv)
{
  unsigned int seed;
  int Ntable, Number;
  double Dratio, Hratio, R, Hmax, DispX, DispZ, fJ, Lx, Ly;
  std::string outfile, config, modfile, modelType;
  bool Mu;

  // Parse command line
  //
  std::string message = "Generate slab initial conditions\n";

  cxxopts::Options options(argv[0], message);

  options.add_options()
    ("h,help",     "Print this help message")
    ("N,number",   "Number of bodies",
     cxxopts::value<int>(Number)->default_value("10000"))
    ("n,ntable",   "Number of points in model table",
     cxxopts::value<int>(Ntable)->default_value("400"))
    ("m,model",    "Model type (LowIso, Sech2, Sech2mu, Sech2Halo)",
     cxxopts::value<std::string>(modelType)->default_value("Sech2"))
    ("d,dratio",    "Ratio of disk to halo density",
     cxxopts::value<double>(Dratio)->default_value("3.0"))
    ("r,hratio",    "Ratio of halo to disk scale height",
     cxxopts::value<double>(Hratio)->default_value("10.0"))
    ("R,vratio",    "Ratio of vertical to horizontal velocity dispersion",
     cxxopts::value<double>(R)->default_value("1.0"))
    ("H,hmax",      "Maximum vertical size in scale heights",
     cxxopts::value<double>(Hmax)->default_value("10.0"))
    ("x,Lx",        "Slab length in the x-direction",
     cxxopts::value<double>(Lx)->default_value("1.0"))
    ("y,Ly",        "Slab length in the y-direction",
     cxxopts::value<double>(Ly)->default_value("1.0"))
    ("X,DispX",     "In-plane velocity variance",
     cxxopts::value<double>(DispX)->default_value("1.0"))
    ("F,fJ",        "Ratio of Jeans length to box scale",
     cxxopts::value<double>(fJ)->default_value("1.0"))
    ("M,Mu",        "Surface density norm for Sech2Halo",
     cxxopts::value<bool>(Mu)->default_value("true"))
    ("c,config",    "Config file",
     cxxopts::value<std::string>(config))
    ("s,seed",      "Random # seed",
     cxxopts::value<unsigned int>(seed))
    ("i,modfile",   "Slab model file for LowIso",
     cxxopts::value<std::string>(modfile)->default_value("slab.model"))
    ("o,outfile",   "output file prefix",
     cxxopts::value<std::string>(outfile)->default_value("slab.bods"))
     ;
  
  cxxopts::ParseResult vm;

  // Parse the command line
  //
  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    return -1;
  }
  
  // Get help
  //
  if (vm.count("help")) {
    std::cout << std::endl << options.help() << std::endl;
    return 0;
  }

  std::ofstream out(outfile);
  if (!out) {
    std::cerr << "Can't open <" << outfile << ">" << std::endl;
    exit(-1);
  }
  out.precision(6);
  out.setf(ios::scientific);

  // Define model
  //
  double h = 1.0;
  DispZ = DispX*R*R;
  std::shared_ptr<OneDModel> model;

  if (modelType.compare("LowIso") == 0) {
    model = std::make_shared<LowIso>(modfile);
  }
  else if (modelType.compare("Sech2mu") == 0) {
    auto tmp = std::make_shared<Sech2mu>(DispZ, DispX);
    h = tmp->get_scale_height();
    if (Hmax>0) tmp->set_hmax(Hmax);
    model = tmp;
  }
  else if (modelType.compare("Sech2") == 0) {
    auto tmp = std::make_shared<Sech2>(DispZ);
    h = tmp->get_scale_height();
    if (Hmax>0) tmp->set_hmax(Hmax);
    model = tmp;
  }
  else if (modelType.compare("Sech2Halo") == 0) {
    Sech2Halo::MU = false;
    if (Mu) Sech2Halo::MU = true;
    auto tmp = std::make_shared<Sech2Halo>(DispZ, Dratio, Hratio);
    h = tmp->get_scale_height();
    if (Hmax>0) tmp->set_hmax(Hmax);
    model = tmp;
  }
  else {
    std::cerr << "non-existent model: " << modelType << std::endl;
    exit(-1);
  }

  // Jeans' length for selecting scale height
  // 
  double maxZ = model->get_max_radius();
  double mu   = model->get_mass(maxZ);
  double KJ   = 2.0*M_PI*mu/DispX;
  double LJ   = 2.0*M_PI/KJ;

  std::cout.setf(ios::left);
  char prev = cout.fill('.');

  std::cout << std::setw(40) << "Model type" << modelType << std::endl;
  std::cout << std::setw(40) << "Surface mass density" << mu << std::endl;
  std::cout << std::setw(40) << "Jeans' length" << LJ << std::endl;
  std::cout << std::setw(40) << "Scale height" << h << std::endl;
  std::cout << std::setw(40) << "Maximum thickness" << maxZ << std::endl;
  
  cout.fill(prev);

  // Make mass table
  //
  std::vector<double> Z(Ntable);
  std::vector<double> M(Ntable);
  double z, dz = 2.0*maxZ/(Ntable-1.0);
  for (int i=0; i<Ntable; i++) {
    Z[i] = -maxZ + dz*i;
    M[i] = model->get_mass(Z[i]);
  }

  std::random_device rd;
  if (vm.count("seed")==0) seed = rd();
  std::mt19937 gen(seed);

  std::uniform_real_distribution<> Unit(0.0, 1.0);
  std::normal_distribution Vv{0.0, sqrt(DispZ)};
  std::normal_distribution Vh{0.0, sqrt(DispX)};

  // Header line
  out << std::setw(10) << Number << std::setw(15) << 0 << std::setw(15) << 0 << std::endl;

  double KE = 0.0;
  double VC = 0.0;
  double mass = mu/Number;

  Eigen::MatrixXd posvel(Number, 6);

#pragma omp parallel for reduction(+:KE,VC)
  for (int n=0; n<Number; n++) {

    posvel.row(n) <<
      Lx*Unit(gen), Ly*Unit(gen), odd2(Unit(gen)*M.back(), M, Z),
      Vh(gen), Vh(gen), Vv(gen);

    KE += posvel(n, 5)*posvel(n, 5);
    VC += posvel(n, 2)*model->get_dpot(posvel(n, 2));
  }
  
  for (int n=0; n<Number; n++) {
    out << std::setw(15) << mass;
    for (int j=0; j<6; j++) out << std::setw(15) << posvel(n, j);
    out << std::endl;
  }

  std::cout << std::endl
	    << "Virial parameters: KE=" << 0.5*mass*KE
	    << " VC=" << mass*VC
	    << " 2T/VC=" << KE/VC << std::endl;

  return 0;
}

