#include <iostream>
#include <iomanip>
#include <strstream>
#include <string>

#include <vector>

#include <math.h>
#include <getopt.h>		// For long options

#include <SLSphere.H>		// Defines biorthogonal SL class
#include <gaussQ.h>		// Gauss-Legendre quadrature

//===========================================================================

				// MPI global variables
				// so one can link to exp libraries
int numprocs, myid, proc_namelen;
char processor_name[MPI_MAX_PROCESSOR_NAME];

MPI_Comm MPI_COMM_SLAVE;

char threading_on = 0;
pthread_mutex_t mem_lock;

//===========================================================================

				// Local function defs
Vector scalar_prod(ScalarType type, double rmin, double rmax, int l, int m,
		   AxiSymBiorth& s, double (*func)(double, int, int), 
		   int numc, int numg);

double plgndr(int l, int m, double costh); 

				// Name of spherical model table
const string model_file_name = "SLGridSph.model";

//===========================================================================

/**
   Class to compute multipole profiles for shift of spherical model
*/
class Shift 
{
private:
  SphericalModelTable *model;
  vector<double> **multi;
  double dR;
  double rmin, rmax;
  int lmax, numr;

public:

  // Constructor
  Shift(double offset, double rmin, double rmax, 
	int lmax, int numr, int numt, int nump);

  // Destructor
  ~Shift();

  // Interpolate on grid
  double shift(double r, int l, int m);

};


Shift::~Shift()
{
  delete model;
}

Shift::Shift(double offset, double Rmin, double Rmax, int Lmax,
	     int numR, int numt, int nump)
{
  rmin = Rmin;
  rmax = Rmax;
  lmax = Lmax;
  numr = numR;

  SphericalModelTable model(model_file_name);
  
  dR = (log(rmax) - log(rmin))/numr;
  double dP = M_PI/nump;
  LegeQuad lq(numt);

  multi = new vector<double>*[lmax+1];
  for (int l=0; l<=lmax; l++) {
    multi[l] = new vector<double>[lmax+1];
    for (int m=0; m<=lmax; m++) 
      multi[l][m] = vector<double>(numr, 0.0);
  }

  double r, theta, phi, rho, facp;
  double x, y, z, R;

  for (int i=0; i<numr; i++) {

    r = rmin*exp(dR*i);
    
    for (int j=0; j<nump; j++) {

      phi = (0.5 + j)*dP;
      
      for (int k=1; k<=numt; k++) {
	
	theta = acos(2.0*(lq.knot(k)-0.5));
	
	x = r * sin(theta)*cos(phi) - offset;
	y = r * sin(theta)*sin(phi);
	z = r * cos(theta);
	
	R = sqrt(x*x + y*y + z*z);

	// 4Pi needed for correct sol'n to Poisson's eqn
	rho = 4.0 * M_PI * model.get_density(R);

	for (int l=0; l<=lmax; l++) {
	  for (int m=0; m<=l; m++) {
	    if (m) facp = 2.0;
	    else facp = 1.0;
	    multi[l][m][i] += facp*dP * 2.0*lq.weight(k) *
	      (0.5*l+0.25)/M_PI*exp(lgamma(1.0+l-m)-lgamma(1.0+l+m)) *
	      plgndr(l, m, 2.0*(lq.knot(k)-0.5)) * cos(phi*m) * rho;
	  }
	}
      }
    }
  }

}


double Shift::shift(double r, int l, int m)
{
  if (abs(m)>l || l>lmax) return 0.0;

  r = max<double>(rmin, r);
  r = min<double>(rmax, r);

  double c[2], ans=0.0;
  int indx;

  indx = (int)( (log(r) - log(rmin))/dR );
  c[1] = (log(r) - log(rmin) - dR*indx);
  c[0] = 1.0 - c[1];
  
  for (int i=0; i<2; i++) ans += c[i]*multi[l][m][indx+i];
  
  return ans;
}


Shift *shift;

double shift_func(double r, int l, int m)
{
  shift->shift(r, l, m);
}

//===========================================================================

/**
   Class to provide reconstructed profiles for given multipole 
   coefficients
*/
class Reconstruct
{
private:
  AxiSymBiorth *biorth;
  Vector **coefs;
  int lmax, nmax;

public:

  // Constructor
  Reconstruct(AxiSymBiorth *bio, double rmin, double rmax, int Lmax, int Nmax);

  // Destructor
  ~Reconstruct();

  // Evaluation of density profile
  double density_eval(double x, double y, double z);

  // Evaluation of potential profile
  double potential_eval(double x, double y, double z);

  // Dump coefficients
  void dump_coefficients(string filename);

};

Reconstruct::~Reconstruct()
{
  for (int l=0; l<=lmax; l++) delete [] coefs[l];
  delete [] coefs;
}

Reconstruct::Reconstruct(AxiSymBiorth *bio, double rmin, double rmax,
			 int Lmax, int Nmax)
{
  biorth = bio;
  lmax = Lmax;
  nmax = Nmax;

  biorth = bio;
  coefs = new Vector*[lmax+1];
  for (int l=0; l<=lmax; l++) {
    coefs[l] = new Vector [lmax+1];
    for (int m=0; m<=l; m++)
      coefs[l][m] = scalar_prod(density, rmin, rmax, l, m, *biorth,
				shift_func, nmax, 200);
  }

}

double Reconstruct::density_eval(double x, double y, double z)
{
  double r = sqrt(x*x + y*y + z*z);
  double costh = z/(r+1.0e-18);
  double phi = atan2(y, x);

  double ans = 0.0;
  for (int l=0; l<=lmax; l++) {
    for (int m=0; m<=l; m++) {
      ans += -biorth->get_dens(r, l, coefs[l][m]) * 
	plgndr(l, m, costh) * 2.0 * cos(phi*m) / (4.0*M_PI);
    }
  }
    
  return ans;
}


double Reconstruct::potential_eval(double x, double y, double z)
{
  double r = sqrt(x*x + y*y + z*z);
  double costh = z/(r+1.0e-18);
  double phi = atan2(y, x);

  double ans = 0.0;
  for (int l=0; l<=lmax; l++) {
    for (int m=0; m<=l; m++) {
      ans += biorth->get_potl(r, l, coefs[l][m]) * 
	plgndr(l, m, costh) * 2.0 * cos(phi*m);
    }
  }
    
  return ans;
}

void Reconstruct::dump_coefficients(string file)
{
  ofstream out(file.c_str());
  if (!out) {
    cerr << "Reconstruct::dump_coefficients: can't open <"
	 << file << "> . . . \n";
    return;
  }
  out.setf(ios::scientific);
  out.precision(8);

  double fac;

  // Header
  out << "# Cosine coefficients only" << endl;
  out << "#" << setw(4) << "l" << setw(5) << "m";
  for (int i=1; i<=nmax; i++) {
    ostrstream ostr;
    ostr << "n=" << i << '\0';
    out << setw(18) <<  ostr.str();
  }
  out << endl;

  // Coefficients
  for (int l=0; l<=lmax; l++) {
    for (int m=0; m<=l; m++) {

      out << setw(5) << l
	  << setw(5) << m;

      fac = sqrt((0.5*l+0.25)/M_PI*
		 exp(lgamma(1.0+l-m)-lgamma(1.0+l+m)));
      if (m) fac *= M_SQRT2;
      
      for (int i=1; i<=nmax; i++) out << setw(18) << coefs[l][m][i]/fac;
      
      out << endl;
    }
  }
}

//===========================================================================

void usage(char *prog)
{
  cout << "Usage:\n\n"
       << prog << " [options]\n\n"
       << setw(15) << "Option" << setw(10) << "Argument" << setw(10) << " " << setw(-40) << "Description" << "\n"
       << "\n"
       << setw(15) << "-m or --mpi" << setw(10) << "No" << setw(10) << " " << setw(-40) << "Turn on MPI for SL computation\n"
       << setw(15) << "-c or --cmap" << setw(10) << "No" << setw(10) << " " << setw(-40) << "Use mapped rather than linear coordinates\n"
       << setw(15) << "--cbio" << setw(10) << "No" << setw(10) << " " << setw(-40) << "Print out basis, if desired\n"
       << setw(15) << "--coefs" << setw(10) << "No" << setw(10) << " " << setw(-40) << "Dump coefficients, if desired\n"
       << setw(15) << "--surface" << setw(10) << "No" << setw(10) << " " << setw(-40) << "Print out surface plots (SM 'ch' style)\n"
       << setw(15) << "--numr" << setw(10) << "Yes" << setw(10) << " " << setw(-40) << "Number of points in radial table\n"
       << setw(15) << "--lmax" << setw(10) << "Yes" << setw(10) << " " << setw(-40) << "Lmax (spherical harmonic expansion)\n"
       << setw(15) << "--nmax" << setw(10) << "Yes" << setw(10) << " " << setw(-40) << "Nmax (radial basis function expansion)\n"
       << setw(15) << "--rmin" << setw(10) << "Yes" << setw(10) << " " << setw(-40) << "Minimum radius for SL basis\n"
       << setw(15) << "--rmax" << setw(10) << "Yes" << setw(10) << " " << setw(-40) << "Maximum radius for SL basis\n"
       << setw(15) << "--rs"<< setw(10) << "Yes" << setw(10) << " " << setw(-40) << "Scale length for radial coordinate mapping\n"
       << setw(15) << "--delr" << setw(10) << "Yes" << setw(10) << " " << setw(-40) << "X-axis offset multipole expansions\n"
       << setw(15) << "--xmax" << setw(10) << "Yes" << setw(10) << " " << setw(-40) << "Length of \"box\" for output profiles\n"
       << setw(15) << "--numx" << setw(10) << "Yes" << setw(10) << " " << setw(-40) << "Number pts for output profiles\n"
       << setw(15) << "--numt" << setw(10) << "Yes" << setw(10) << " " << setw(-40) << "Number knots for cos(theta) integral\n"
       << setw(15) << "--nump" << setw(10) << "Yes" << setw(10) << " " << setw(-40) << "Number knots for phi integral\n"
       << setw(15) << "--file" << setw(10) << "Yes" << setw(10) << " " << setw(-40) << "File name prefix for output\n"
       << "\n";

  exit(0);
}

int 
main(int argc, char** argv)
{
				// Default values defined here
  bool use_mpi = false;
  bool check_bio = false;
  bool surface = false;
  bool dump = false;
  int cmap = 0;
  double scale = 1.0;

  int Lmax=2, nmax=10;
  int numr=1000;
  double rmin=0.001, rmax=1.95, rs=0.067;
  double delr=0.01, xmax=1.0;
  int numx=100, numt = 20, nump = 20;

  string outfile = "slshift";

  int c;
  while (1) {
    int this_option_optind = optind ? optind : 1;
    int option_index = 0;
    static struct option long_options[] = {
      {"mpi", 0, 0, 0},		// Turn on MPI for SL computation
      {"cmap", 0, 0, 0},	// Use mapped rather than linear coordinates
      {"cbio", 0, 0, 0},	// Print out basis, if desired
      {"coefs", 0, 0, 0},	// Dump coefficients, if desired
      {"surface", 0, 0, 0},	// Print out surface plots (SM 'ch' style)
      {"numr", 1, 0, 0},	// Number of points in radial table
      {"lmax", 1, 0, 0},	// Lmax (spherical harmonic expansion)
      {"nmax", 1, 0, 0},	// Nmax (radial basis function expansion)
      {"rmin", 1, 0, 0},	// Minimum radius for SL basis
      {"rmax", 1, 0, 0},	// Maximum radius for SL basis
      {"rs", 1, 0, 0},		// Scale length for radial coordinate mapping
      {"delr", 1, 0, 0},	// X-axis offset multipole expansions
      {"xmax", 1, 0, 0},	// Length of "box" for output profiles
      {"numx", 1, 0, 0},	// Number pts for output profiles
      {"numt", 1, 0, 0},	// Number knots for cos(theta) integral
      {"nump", 1, 0, 0},	// Number knots for phi integral
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
	} else if (!optname.compare("cbio")) {
	  check_bio = true;
	} else if (!optname.compare("cmap")) {
	  cmap = 1;
	} else if (!optname.compare("surface")) {
	  surface = true;
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
	} else if (!optname.compare("nump")) {
	  nump = atoi(optarg);
	} else if (!optname.compare("numt")) {
	  numt = atoi(optarg);
	} else if (!optname.compare("file")) {
	  outfile = string(optarg);
	} else {
	  cout << "Option " << long_options[option_index].name;
	  if (optarg) cout << " with arg " << optarg;
	  cout << " is not defined \n";
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
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &proc_namelen);
  }

  if (use_mpi) {

    SLGridSph::mpi = 1;		// Turn on MPI

  } else {

    SLGridSph::mpi = 0;		// Turn off MPI
  }

  //===================
  // Set up basis
  //===================

				// Generate Sturm-Liouville grid
  SLSphere ortho(Lmax, nmax, numr, rmin, rmax, cmap, rs);
  
  if (check_bio) {
    ofstream out("slshift.bio");

    cout << "Number of points? ";
    int num;
    cin >> num;

    cout << "L, N? ";
    int L, N;
    cin >> L;
    cin >> N;
    
    double ximin = ortho.r_to_rb(rmin);
    double ximax = ortho.r_to_rb(rmax);

    double x, r;

    for (int i=0; i<num; i++) {
      x = ximin + (ximax-ximin)*(0.5 + i)/num;
      r = ortho.rb_to_r(x);
      out << setw(15) << x
	  << setw(15) << r
	  << setw(15) << ortho.potl(N, L, x)
	  << setw(15) << ortho.dens(N, L, x)
	  << setw(15) << ortho.potlR(N, L, r)
	  << endl;
    }
  }

  if (!use_mpi || myid==0) {

    //===================
    // Get coefficients
    //===================

				// Construct coefficients
    shift = new Shift(delr, rmin, rmax, Lmax, numr, numt, nump);

    //===================
    // Shifted profile
    //===================
				// Reconstruction
    Reconstruct recon(&ortho, rmin, rmax, Lmax, nmax);

    string ostr(outfile);
    ostr += ".profile";
    ofstream out(ostr.c_str());

    double x, dx = 2.0*xmax/(numx-1);
    for (int i=0; i<numx; i++) {
      x = -xmax + dx*i;
      cout << setw(15) << x
	   << setw(15) << recon.density_eval(x, 0.0, 0.0)
	   << setw(15) << recon.potential_eval(x, 0.0, 0.0)
	   << setw(15) << recon.density_eval(delr, x, 0.0)
	   << setw(15) << recon.potential_eval(delr, x, 0.0)
	   << endl;
      out  << setw(15) << x
	   << setw(15) << recon.density_eval(x, 0.0, 0.0)
	   << setw(15) << recon.potential_eval(x, 0.0, 0.0)
	   << setw(15) << recon.density_eval(delr, x, 0.0)
	   << setw(15) << recon.potential_eval(delr, x, 0.0)
	   << endl;
    }

    if (surface) {
      double x, y, dxy = 2.0*xmax/(numx-1);
      float z;

      ofstream *out = new ofstream [2];
      string suffix[2] = {".potl\0", ".dens\0"};
      for (int i=0; i<2; i++) {
	string ostr(outfile);
	ostr += suffix[i];
	out[i].open(ostr.c_str());
	for (int j=0; j<2; j++) out[i].write(&numx, sizeof(int));
	for (int j=0; j<2; j++) {
	  out[i].write(&(z=-xmax), sizeof(float));
	  out[i].write(&(z= xmax), sizeof(float));
	}
      }
      
      for (int j=0; j<numx; j++) {
	y = -xmax + dxy*j;
	for (int i=0; i<numx; i++) {
	  x = -xmax + dxy*i;
	  
	  out[0].write(&(z=recon.potential_eval(x, y, 0.0)), sizeof(float));
	  out[1].write(&(z=recon.density_eval(x, y, 0.0)), sizeof(float));
	}
      }
    }

    if (dump) {
      string ostr(outfile);
      ostr += ".coefs";
      recon.dump_coefficients(ostr);
    }

  }

  return 0;
}
