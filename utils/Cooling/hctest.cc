
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <random>
#include <cmath>

using namespace std;

#include <Timer.H>
#include <HeatCool.H>
#include <cxxopts.H>		// Option parsing
#include <libvars.H>		// EXP global variables

int main(int argc, char** argv)
{
  double Nmin=1.0e-8, Nmax=1.0e10, Tmin=1000, Tmax=1e7;
  unsigned Nnum=100, Tnum=100, nnum=40, tnum=40;
  string filename = "hctest.dat";
  string cachefile = ".HeatCool";

  cxxopts::Options options("hctest", "Heating and cooling test");

  options.add_options()
    ("h,help", "print this help message")
    ("1,Nmin", "minimum number density (H/cc)",
     cxxopts::value<double>(Nmin)->default_value("1.0e-08"))
    ("2,Nmax", "maximum number density (H/cc)",
     cxxopts::value<double>(Nmax)->default_value("1.0e10"))
    ("3,Tmin", "minimum temperature",
     cxxopts::value<double>(Tmin)->default_value("1.0e+02"))
    ("4,Tmax", "maximum temperature",
     cxxopts::value<double>(Tmax)->default_value("1.0e+07"))
    ("f,filename", "output filename",
     cxxopts::value<std::string>(filename)->default_value("hctest.dat"))
    ("c,cachefile", "cache file",
     cxxopts::value<std::string>(cachefile)->default_value(".HeatCool"))
    ("Nnum", "number of density points",
     cxxopts::value<unsigned int>(Nnum)->default_value("100"))
    ("Tnum", "number of temperature points",
     cxxopts::value<unsigned int>(Tnum)->default_value("100"))
    ("nnum", "number of density points for grid evaluation",
     cxxopts::value<unsigned int>(nnum)->default_value("40"))
    ("Tnum", "number of temperature points grid evaluation",
     cxxopts::value<unsigned int>(tnum)->default_value("40"))
    ;

  

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    std::cout << "Option error: " << e.what() << std::endl;
    exit(-1);
  }

  if (vm.count("help")) {
    std::cout << options.help() << std::endl;
    exit(-1);
  }

  //=====================
  // Setup grid
  //=====================

  HeatCool hc(Nmin, Nmax, Tmin, Tmax, Nnum, Tnum, cachefile);

  //=====================
  // Output in GP format
  //=====================

  ofstream out(filename.c_str());
  if (!out) {
    cerr << "Could not open <" << filename << "> for output"
	 << endl;
    return -1;
  }

  double dN = (log(Nmax) - log(Nmin))/(nnum-1);
  double dT = (log(Tmax) - log(Tmin))/(tnum-1);

  double N, T;

  for (unsigned n=0; n<nnum; n++) {

    N = Nmin*exp(dN*n);

    for (unsigned t=0; t<tnum; t++) {

      T = Tmin*exp(dT*t);
      hc.setPoint(N, T);
      out << setw(14) << N
	  << setw(14) << T
	  << setw(14) << hc.CoolRate()
	  << endl;
    }
    out << endl;
  }

  out.close();

  //=====================
  // Timing
  //=====================
  Timer one, two;
  double oneSoFar, twoSoFar;
  double tst;

  one.start();
  for (unsigned n=0; n<nnum; n++) {
    N = Nmin*exp(dN*n);
    for (unsigned t=0; t<tnum; t++) {
      T = Tmin*exp(dT*t);
      tst = hc.CoolRate(N, T);
    }
  }
  oneSoFar = one.stop();

  two.start();
  for (unsigned n=0; n<nnum; n++) {
    N = Nmin*exp(dN*n);
    for (unsigned t=0; t<tnum; t++) {
      T = Tmin*exp(dT*t);
      HeatCool tmp(N, T);
    }
  }
  twoSoFar = two.stop();

  cout << "Interpolate = " << oneSoFar << endl;
  cout << "Computation = " << twoSoFar << endl;

  //=====================
  // Test accuracy
  //=====================

  double interp, exact, maxrel=0.0, maxabs=0.0;
  double worst_rel[4], worst_abs[4];
  for (unsigned n=0; n<nnum; n++) {
    N = Nmin*exp(dN*n);
    for (unsigned t=0; t<tnum; t++) {
      T = Tmin*exp(dT*t);
      interp = hc.CoolRate(N, T);
      HeatCool tmp(N, T);
      exact = tmp.CoolRate();

      tst = fabs(interp-exact);
      if (maxabs < tst) {
	maxabs = tst;
	worst_abs[0] = N;
	worst_abs[1] = T;
	worst_abs[2] = exact;
	worst_abs[3] = interp;
      }

      if (exact>0.0) {
	tst /= exact;
	if (maxrel < tst) {
	  maxrel = tst;
	  worst_rel[0] = N;
	  worst_rel[1] = T;
	  worst_rel[2] = exact;
	  worst_rel[3] = interp;
	}
      }
    }
  }

  cout << "Maximum relative error = " << maxrel << endl;
  cout << "    " 
       << setw(14) << worst_rel[0]
       << setw(14) << worst_rel[1]
       << setw(14) << worst_rel[2]
       << setw(14) << worst_rel[3]
       << endl;
  cout << "Maximum absolute error = " << maxabs << endl;
  cout << "    " 
       << setw(14) << worst_abs[0]
       << setw(14) << worst_abs[1]
       << setw(14) << worst_abs[2]
       << setw(14) << worst_abs[3]
       << endl;

  return 0;
}




