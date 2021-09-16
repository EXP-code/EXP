#include <math.h>
#include <getopt.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#include <boost/random/mersenne_twister.hpp>

using namespace std;

#include <Timer.H>
#include <HeatCool.H>

int myid = 0;
char threading_on = 0;
string outdir = "";
string runtag = "test";
pthread_mutex_t mem_lock;
boost::mt19937 random_gen;

//===========================================================================

void usage(char *prog)
{
  static unsigned f1=20, f2=12, f3=40;

  cout << "Usage:" << endl << endl
       << prog << " [options]" << endl << endl
       << setiosflags(ios::left)
       << setw(f1) << "Option" << setw(f2) << "Argument" << setw(10) << " " 
       << setw(f3) << "Description" << endl
       << endl
       << setw(f1) << "-1 or --Nmin" << setw(f2) << "0.000001" << setw(10) << " " 
       << setw(f3) << "Minimum number density (H/cc)" << endl
       << setw(f1) << "-c or --cache" << setw(f2) << ".HeatCool" << setw(10) << " " 
       << setw(f3) << "HeatCool cache filename" << endl
       << setw(f1) << "-f or --file" << setw(f2) << "hctest.dat" << setw(10) << " " 
       << setw(f3) << "Output data" << endl
       << resetiosflags(ios::left)
       << endl;

  exit(0);
}

int main(int argc, char** argv)
{
  double Nmin=1.0e-8, Nmax=1.0e10, Tmin=1000, Tmax=1e7;
  unsigned Nnum=100, Tnum=100, nnum=40, tnum=40;
  string filename = "hctest.dat";
  string cachefile = ".HeatCool";

  int c;
  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
      {"Nmin", 1, 0, 0},
      {"Nmax", 1, 0, 0},
      {"Tmin", 1, 0, 0},
      {"Tmax", 1, 0, 0},
      {"Nnum", 1, 0, 0},
      {"Tnum", 1, 0, 0},
      {"nnum", 1, 0, 0},
      {"tnum", 1, 0, 0},
      {"cache", 1, 0, 0},
      {"file", 1, 0, 0},
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "1:2:3:4:N:T:n:t:c:f:h",
		     long_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 0:
      {
	string optname(long_options[option_index].name);

	if (!optname.compare("Nmin")) {
	  Nmin = atof(optarg);
	} else if (!optname.compare("Nmax")) {
	  Nmax = atoi(optarg);
	} else if (!optname.compare("Tmin")) {
	  Tmin = atoi(optarg);
	} else if (!optname.compare("Tmax")) {
	  Tmax = atoi(optarg);
	} else if (!optname.compare("Nnum")) {
	  Nnum = atoi(optarg);
	} else if (!optname.compare("Tnum")) {
	  Tnum = atoi(optarg);
	} else if (!optname.compare("nnum")) {
	  nnum = atoi(optarg);
	} else if (!optname.compare("tnum")) {
	  tnum = atoi(optarg);
	} else if (!optname.compare("cache")) {
	  cachefile = optarg;
	} else if (!optname.compare("file")) {
	  filename = optarg;
	} else {
	  cout << "Option " << long_options[option_index].name;
	  if (optarg) cout << " with arg " << optarg;
	  cout << " is not defined " << endl;
	  exit(0);
	}
      }
      break;

    case '1':
      Nmin = atof(optarg);
      break;

    case '2':
      Nmax = atof(optarg);
      break;

    case '3':
      Tmin = atof(optarg);
      break;

    case '4':
      Tmax = atof(optarg);
      break;

    case 'N':
      Nnum = atoi(optarg);
      break;

    case 'T':
      Tnum = atoi(optarg);
      break;

    case 'n':
      nnum = atoi(optarg);
      break;

    case 't':
      tnum = atoi(optarg);
      break;

    case 'f':
      filename = optarg;
      break;

    case 'c':
      cachefile = optarg;
      break;

    case 'h':
    default:
      usage(argv[0]);
    }

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




