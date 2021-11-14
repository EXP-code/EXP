#include <iostream>
#include <iomanip>
#include <sstream>
#include <memory>
#include <string>
#include <vector>
#include <cmath>

#include <getopt.h>		// For long options

#include <global.H>		// EXP globals
#include <localmpi.H>		// MPI globals
#include <SLSphere.H>		// Defines biorthogonal SL class
#include <CylindricalDisk.H>	// The axisymmetric potential solver
#include <gaussQ.H>		// Gauss-Legendre quadrature
#include <exponential.H>	// Exponential disk
#include <toomre.H>		// Toomre disk (including Kuzmin)

//===========================================================================

void usage(char *prog)
{
  cout << "Usage:" << endl << endl
       << prog << " [options]" << endl << endl
       << setw(15) << "Option" << setw(10) << "Argument" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Description" << endl << resetiosflags(ios::left)
       << endl
       << setw(15) << "-m or --mpi" << setw(10) << "No" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Turn on MPI for SL computation" << endl << resetiosflags(ios::left)
       << setw(15) << "-c or --cmap" << setw(10) << "No" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Use mapped rather than linear coordinates" << endl << resetiosflags(ios::left)
       << setw(15) << "--Kuzmin" << setw(10) << "No" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Use the Kuzmin disk surface density" << endl << resetiosflags(ios::left)
       << setw(15) << "--coefs" << setw(10) << "No" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Dump coefficients, if desired" << endl << resetiosflags(ios::left)
       << setw(15) << "--numr" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Number of points in radial table" << endl << resetiosflags(ios::left)
       << setw(15) << "--lmax" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Lmax (spherical harmonic expansion)" << endl << resetiosflags(ios::left)
       << setw(15) << "--nmax" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Nmax (radial basis function expansion)" << endl << resetiosflags(ios::left)
       << setw(15) << "--rmin" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Minimum radius for SL basis" << endl << resetiosflags(ios::left)
       << setw(15) << "--rmax" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Maximum radius for SL basis" << endl << resetiosflags(ios::left)
       << setw(15) << "--delr" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "X-axis offset multipole expansions" << endl << resetiosflags(ios::left)
       << setw(15) << "--delta" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Fractional offset for difference derivatives" << endl << resetiosflags(ios::left)
       << setw(15) << "--xmax" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Length of \"box\" for output profiles" << endl << resetiosflags(ios::left)
       << setw(15) << "--zmax" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Height of \"box\" for output profiles" << endl << resetiosflags(ios::left)
       << setw(15) << "--numx" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Number pts in length for output profiles" << endl << resetiosflags(ios::left)
       << setw(15) << "--numz" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Number pts in height for output profiles" << endl << resetiosflags(ios::left)
       << setw(15) << "--numt" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Number knots for cos(theta) integral" << endl << resetiosflags(ios::left)
       << setw(15) << "--numg" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "Number of points for interpolation grid" << endl << resetiosflags(ios::left)
       << setw(15) << "--file" << setw(10) << "Yes" << setw(10) << " " << setiosflags(ios::left) << setw(40) << "File name prefix for output" << endl << resetiosflags(ios::left)
       << endl;

  exit(0);
}

string outdir, runtag;
std::mt19937 random_gen;

int 
main(int argc, char** argv)
{
				// Default values defined here
  bool use_mpi = false;
  bool Kuzmin  = false;
  bool dump    = false;
  int cmap     = 0;
  double scale = 1.0;

  int Lmax=16, Nmax=10;
  int numr=10000;
  double rmin=0.0001, rmax=1.0;
  double delr=0.01, xmax=1.0, zmax=0.1;
  int numx=100, numz=100, numt=1000, numg= 100;

  string outfile = "diskpot";

  int c;
  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
      {"mpi", 0, 0, 0},		// Turn on MPI for SL computation
      {"cmap", 0, 0, 0},	// Use mapped rather than linear coordinates
      {"Kuzmin", 0, 0, 0},	// Use a Kuzmin rather than Exponential disk
      {"coefs", 0, 0, 0},	// Dump coefficients, if desired
      {"numr", 1, 0, 0},	// Number of points in radial table
      {"lmax", 1, 0, 0},	// Lmax (spherical harmonic expansion)
      {"nmax", 1, 0, 0},	// Nmax (radial basis function expansion)
      {"rmin", 1, 0, 0},	// Minimum radius for SL basis
      {"rmax", 1, 0, 0},	// Maximum radius for SL basis
      {"delr", 1, 0, 0},	// X-axis offset multipole expansions
      {"delta", 1, 0, 0},	// Fractional offset for difference derivs
      {"xmax", 1, 0, 0},	// Length of "box" for output profiles
      {"numx", 1, 0, 0},	// Number pts in width for output profiles
      {"zmax", 1, 0, 0},	// Height of "box" for output profiles
      {"numz", 1, 0, 0},	// Number pts in height for output profiles
      {"numt", 1, 0, 0},	// Number knots for cos(theta) integral
      {"numg", 1, 0, 0},	// Number points for interpolation grid
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
	} else if (!optname.compare("Kuzmin")) {
	  Kuzmin = true;
	} else if (!optname.compare("coefs")) {
	  dump = true;
	} else if (!optname.compare("numr")) {
	  numr = atoi(optarg);
	} else if (!optname.compare("lmax")) {
	  Lmax = atoi(optarg);
	} else if (!optname.compare("nmax")) {
	  Nmax = atoi(optarg);
	} else if (!optname.compare("rmin")) {
	  rmin = atof(optarg);
	} else if (!optname.compare("rmax")) {
	  rmax = atof(optarg);
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
	} else if (!optname.compare("numt")) {
	  numt = atoi(optarg);
	} else if (!optname.compare("numg")) {
	  numg = atoi(optarg);
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

  vector<double> param(3);
  std::shared_ptr<CylindricalDisk> disk;

  if (Kuzmin) {
    param[0] = 1.0;		// Velocity scale
    param[1] = 0.01;		// Scale length
    param[2] = 0.001;		// Scale height
    
    disk = std::make_shared<KuzminDisk>();

  } else {
    param[0] = 0.1;		// Disk mass
    param[1] = 0.01;		// Scale length
    param[2] = 0.001;		// Scale height

    disk = std::make_shared<CylindricalDisk>();
  }

  disk->Initialize(rmin, rmax, true, Nmax, Lmax, numr, numt, numg, param);

  //=====================
  // Print the expansion
  //=====================

  double x, dx = 2.0*xmax/(numx-1);
  double z, dz = 2.0*zmax/(numz-1);
  
  const int nfiles = 4;
  std::vector<ofstream> out(nfiles);
  string suffix[nfiles] = {".potl", ".dens", ".force", ".force0"};
  for (int i=0; i<nfiles; i++) {
    string ostr(outfile);
    ostr += suffix[i];
    out[i].open(ostr.c_str());
  }
      
  double fx, fz, fx1, fz1, mf;
  double dfx, dfz;
  const double frac = 0.03;

  for (int j=0; j<numz; j++) {
    z = -zmax + dz*j;

    for (int i=0; i<numx; i++) {
      x = -xmax + dx*i;
      
      for (int k=0; k<nfiles-1; k++) out[k] << setw(16) << x << setw(16) << z;
      out[0] << setw(16) << disk->potential_eval(x, 0.0, z);
      out[1] << setw(16) << disk->density_eval(x, 0.0, z);

      //
      // Set the derivative mesh spacing
      //
      dfx = max<double>(1.0e-5, x*frac);
      dfz = max<double>(1.0e-6, z*frac);

      //
      // 3-pt formula
      //
      /*
      fx = (disk->potential_eval(x-dfx, 0.0, z    ) - 
	    disk->potential_eval(x+dfx, 0.0, z    ) )/(2.0*dfx);

      fz = (disk->potential_eval(x    , 0.0, z-dfz) - 
	    disk->potential_eval(x    , 0.0, z+dfz) )/(2.0*dfz);
      */

      //
      // 5-pt formula
      //
      fx = -(
	     -disk->potential_eval(x+2.0*dfx, 0.0, z    )
	     +8.0*disk->potential_eval(x+dfx, 0.0, z    )
	     -8.0*disk->potential_eval(x-dfx, 0.0, z    )
	     +disk->potential_eval(x-2.0*dfx, 0.0, z    )
	     )/(12.0*dfx);

      fz = -(
	     -disk->potential_eval(x    , 0.0, z+2.0*dfz)
	     +8.0*disk->potential_eval(x, 0.0, z+dfz    )
	     -8.0*disk->potential_eval(x, 0.0, z-dfz    )
	     +disk->potential_eval(x    , 0.0, z-2.0*dfz)
	     )/(12.0*dfz);

      disk->force_eval(x, 0.0, z, fx1, fz1);
      if (x<0.0) fx1 *= -1.0;

      out[2] << setw(16) << fx << setw(16) << fz
	     << setw(16) << fx - fx1 << setw(16) << fz - fz1;

      for (int k=0; k<nfiles-1; k++) out[k] << endl;
    }
    for (int k=0; k<nfiles-1; k++) out[k] << endl;
  }
  
  
  std::shared_ptr<AxiSymModel> edisk;
  double A, mass;

  if (Kuzmin) {
    edisk = std::make_shared<ToomreDisk>();
    A = param[1];
    mass = param[0]*param[0];
  } else {
    edisk = std::make_shared<ExponentialDisk>(param[1]);
    A = param[1];
    mass = param[0];
  }

  for (int i=0; i<numx; i++) {

    x = -xmax + dx*i;
      
    //
    // Set the derivative mesh spacing
    //
    dfx = max<double>(1.0e-5, x*frac);

    fx = -(
	   -disk->potential_eval(x+2.0*dfx, 0.0, 0.0)
	   +8.0*disk->potential_eval(x+dfx, 0.0, 0.0)
	   -8.0*disk->potential_eval(x-dfx, 0.0, 0.0)
	   +disk->potential_eval(x-2.0*dfx, 0.0, 0.0)
	   )/(12.0*dfx);

    disk->force_eval(x, 0.0, 0.0, fx1, fz1);
    if (x<0.0) fx1 *= -1.0;

    // Force from 2-d model
    double xx = fabs(x);
    double y = xx/A;
    
    out[3] << setw(16) << x << setw(16) << fx
	   << setw(16) << fx - fx1;
    
    if (Kuzmin) {
      out[3] << setw(16) << -mass/A*edisk->get_dpot(y)*y/x;
    } else {
      out[3] << setw(16) << -mass*edisk->get_dpot(xx)*x/(xx+1.0e-10);
    }
    out[3] << endl;

  }
  
  if (dump) {
    string ostr(outfile);
    ostr += ".coefs";
    disk->dump_coefficients(ostr);
  }
  
  return 0;
}
