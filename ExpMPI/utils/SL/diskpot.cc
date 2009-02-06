#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include <getopt.h>		// For long options

#include <localmpi.h>
#include <SLSphere.H>		// Defines biorthogonal SL class
#include <CylindricalDisk.H>	// The axisymmetric potential solver
#include <gaussQ.h>		// Gauss-Legendre quadrature

//===========================================================================

				// so one can link to exp libraries
char threading_on = 0;
pthread_mutex_t mem_lock;

//===========================================================================


//===========================================================================

void usage(char *prog)
{
  cout << "Usage:" << endl << endl
       << prog << " [options]" << endl << endl
       << setw(15) << "Option" << setw(10) << "Argument" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Description" << endl << resetiosflags(ios::left)
       << endl
       << setw(15) << "-m or --mpi" << setw(10) << "No" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Turn on MPI for SL computation" << endl << resetiosflags(ios::left)
       << setw(15) << "-c or --cmap" << setw(10) << "No" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Use mapped rather than linear coordinates" << endl << resetiosflags(ios::left)
       << setw(15) << "--coefs" << setw(10) << "No" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Dump coefficients, if desired" << endl << resetiosflags(ios::left)
       << setw(15) << "--numr" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Number of points in radial table" << endl << resetiosflags(ios::left)
       << setw(15) << "--lmax" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Lmax (spherical harmonic expansion)" << endl << resetiosflags(ios::left)
       << setw(15) << "--nmax" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Nmax (radial basis function expansion)" << endl << resetiosflags(ios::left)
       << setw(15) << "--rmin" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Minimum radius for SL basis" << endl << resetiosflags(ios::left)
       << setw(15) << "--rmax" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Maximum radius for SL basis" << endl << resetiosflags(ios::left)
       << setw(15) << "--rs"<< setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Scale length for radial coordinate mapping" << endl << resetiosflags(ios::left)
       << setw(15) << "--delr" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "X-axis offset multipole expansions" << endl << resetiosflags(ios::left)
       << setw(15) << "--delta" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Fractional offset for difference derivatives" << endl << resetiosflags(ios::left)
       << setw(15) << "--xmax" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Length of \"box\" for output profiles" << endl << resetiosflags(ios::left)
       << setw(15) << "--zmax" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Height of \"box\" for output profiles" << endl << resetiosflags(ios::left)
       << setw(15) << "--numx" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Number pts in length for output profiles" << endl << resetiosflags(ios::left)
       << setw(15) << "--numz" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Number pts in height for output profiles" << endl << resetiosflags(ios::left)
       << setw(15) << "--numt" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Number knots for cos(theta) integral" << endl << resetiosflags(ios::left)
       << setw(15) << "--file" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "File name prefix for output" << endl << resetiosflags(ios::left)
       << endl;

  exit(0);
}

int 
main(int argc, char** argv)
{
				// Default values defined here
  bool use_mpi = false;
  bool dump = false;
  int cmap = 0;
  double scale = 1.0;

  int Lmax=2, nmax=10;
  int numr=10000;
  double rmin=0.0001, rmax=1.95, rs=0.067;
  double delr=0.01, xmax=5.0, zmax=0.5;
  int numx=100, numz=100, numt = 40, nump = 40;

  string outfile = "slshift";

  int c;
  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
      {"mpi", 0, 0, 0},		// Turn on MPI for SL computation
      {"cmap", 0, 0, 0},	// Use mapped rather than linear coordinates
      {"coefs", 0, 0, 0},	// Dump coefficients, if desired
      {"numr", 1, 0, 0},	// Number of points in radial table
      {"lmax", 1, 0, 0},	// Lmax (spherical harmonic expansion)
      {"nmax", 1, 0, 0},	// Nmax (radial basis function expansion)
      {"rmin", 1, 0, 0},	// Minimum radius for SL basis
      {"rmax", 1, 0, 0},	// Maximum radius for SL basis
      {"rs", 1, 0, 0},		// Scale length for radial coordinate mapping
      {"delr", 1, 0, 0},	// X-axis offset multipole expansions
      {"delta", 1, 0, 0},	// Fractional offset for difference derivs
      {"xmax", 1, 0, 0},	// Length of "box" for output profiles
      {"numx", 1, 0, 0},	// Number pts in width for output profiles
      {"zmax", 1, 0, 0},	// Height of "box" for output profiles
      {"numz", 1, 0, 0},	// Number pts in height for output profiles
      {"numt", 1, 0, 0},	// Number knots for cos(theta) integral
      {"file", 1, 0, 0},	// File name prefix for output
      {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "cmh",
		     long_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 0:
      {
	string optname(long_options[option_index].name);

	if (!optname.compare("mpi")) {
	  use_mpi = true;
	} else if (!optname.compare("cmap")) {
	  cmap = 1;
	} else if (!optname.compare("coefs")) {
	  dump = true;
	} else if (!optname.compare("numr")) {
	  numr = atoi(optarg);
	} else if (!optname.compare("lmax")) {
	  Lmax = atoi(optarg);
	} else if (!optname.compare("nmax")) {
	  nmax = atoi(optarg);
	} else if (!optname.compare("rmin")) {
	  rmin = atof(optarg);
	} else if (!optname.compare("rmax")) {
	  rmax = atof(optarg);
	} else if (!optname.compare("rs")) {
	  rs = atof(optarg);
	} else if (!optname.compare("delr")) {
	  delr = atof(optarg);
	} else if (!optname.compare("xmax")) {
	  xmax = atof(optarg);
	} else if (!optname.compare("numx")) {
	  numx = atoi(optarg);
	} else if (!optname.compare("zmax")) {
	  zmax = atof(optarg);
	} else if (!optname.compare("numz")) {
	  numz = atoi(optarg);
	} else if (!optname.compare("nump")) {
	  nump = atoi(optarg);
	} else if (!optname.compare("numt")) {
	  numt = atoi(optarg);
	} else if (!optname.compare("file")) {
	  outfile = string(optarg);
	} else {
	  cout << "Option " << long_options[option_index].name;
	  if (optarg) cout << " with arg " << optarg;
	  cout << " is not defined" << endl;
	  exit(0);
	}
      }
      break;

    case 'c':
      cmap = 1;
      break;

    case 'm':
      use_mpi = true;
      break;

    case 'h':
    default:
      usage(argv[0]);
    }

  }

  //===================
  // MPI preliminaries 
  //===================
  if (use_mpi) {
    local_init_mpi(argc, argv);
  }

  if (use_mpi) {

    SLGridSph::mpi = 1;		// Turn on MPI

  } else {

    SLGridSph::mpi = 0;		// Turn off MPI
  }

  //===================
  // Set up the disk
  //===================

  CylindricalDisk disk(rmin, rmax, nmax, Lmax, numr, numt);

  //=====================
  // Print the expansion
  //=====================

  double x, dx = 2.0*xmax/(numx-1);
  double z, dz = 2.0*zmax/(numz-1);
  
  ofstream *out = new ofstream [2];
  string suffix[2] = {".potl\0", ".dens\0"};
  for (int i=0; i<2; i++) {
    string ostr(outfile);
    ostr += suffix[i];
    out[i].open(ostr.c_str());
  }
      
  for (int j=0; j<numz; j++) {
    z = -zmax + dz*j;
    for (int i=0; i<numx; i++) {
      x = -xmax + dx*i;
      
      out[0] << setw(16) << x << setw(16) << z 
	     << setw(16) << disk.potential_eval(x, 0.0, z) << endl;
      
      out[1] << setw(16) << x << setw(16) << z 
	     << setw(16) << disk.density_eval(x, 0.0, z) << endl;
    }
  }
  
  if (dump) {
    string ostr(outfile);
    ostr += ".coefs";
    disk.dump_coefficients(ostr);
  }
  
  return 0;
}
