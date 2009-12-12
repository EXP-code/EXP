#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

#include <getopt.h>		// For long options

#include <localmpi.h>
#include <ProgramParam.H>
#include <SLSphere.H>		// Defines biorthogonal SL class
#include <gaussQ.h>		// Gauss-Legendre quadrature


//===========================================================================

bool     use_mpi;
bool     check_bio;
bool     surface;
bool     dump;
int      cmap;
int      Lmax;
int      nmax;
int      numr;
double   rmin;
double   rmax;
double   rs;
double   scale;
double   delr;
double   delta;
double   xmax;
int      numx;
int      numt;
int      nump;
string   outfile;


program_option init[] = {
  {"use_mpi",	"bool",		"false",	"using parallel computation"},
  {"check_bio",	"bool",		"false",	"check consistency of biorthogonal set"},
  {"surface",	"bool",		"false",	"print out surface cuts of the reconstructed shift"},
  {"dump",	"bool",		"false",	"print the coefficients to a file"},
  {"cmap",	"int",		"0",		"coordinate scaling in SphereSL"},
  {"scale",	"double",	"1.0",		"scaling from real coordinates to table"},
  {"Lmax",	"int",		"2",		"maximum number of angular harmonics in the expansion"},
  {"nmax",	"int",		"10",		"maximum number of radial harmonics in the expansion"},
  {"numr",	"int",		"1000",		"radial knots in the shift operator"},
  {"rmin",	"double",	"0.0001",	"minimum radius for the shift operator"},
  {"rmax",	"double",	"1.95",		"maximum radius for the shift operator"},
  {"rs",	"double",	"0.067",	"cmap scale factor"},
  {"delr",	"double",	"0.003",	"horizontal shift for test"},
  {"xmax",	"double",	"1.0",		"radial scale for output"},
  {"numx",	"int",		"100",		"number of output knots per side"},
  {"numt",	"int",		"40",		"number of theta knots in shift operator"},
  {"nump",	"int",		"40",		"number of phi knots in shift operator"},
  {"delta",	"double",	"0.05",		"fractional displacement for two-point numerical derivative"},
  {"outfile",	"string",	"slshift",	"output file prefix"},
  {"\0",	"\0",		"\0",		"\0"}
};


const char *desc = "Check the consistency of a linear shift in the basis";

ProgramParam config(desc, init);


void assign_values()
{
  use_mpi   = config.get<bool>  ("use_mpi");
  check_bio = config.get<bool>  ("check_bio");
  surface   = config.get<bool>  ("surface");
  dump      = config.get<bool>  ("dump");
  cmap      = config.get<int>   ("cmap");
  scale     = config.get<double>("scale");
  Lmax      = config.get<int>   ("Lmax");
  nmax      = config.get<int>   ("nmax");
  numr      = config.get<int>   ("numr");
  rmin      = config.get<double>("rmin");
  rmax      = config.get<double>("rmax");
  rs        = config.get<double>("rs");
  delr      = config.get<double>("delr");
  xmax      = config.get<double>("xmax");
  numx      = config.get<int>   ("numx");
  numt      = config.get<int>   ("numt");
  nump      = config.get<int>   ("nump");
  delta     = config.get<double>("delta");
  outfile   = config.get<string>("outfile");
};


//===========================================================================

				// so one can link to exp libraries
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
	    else   facp = 1.0;
	    multi[l][m][i] += facp*2.0*dP * 2.0*lq.weight(k) *
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
				shift_func, nmax, 400);
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
	plgndr(l, m, costh) * cos(phi*m) / (4.0*M_PI);
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
	plgndr(l, m, costh) * cos(phi*m);
    }
  }
    
  return ans;
}

void Reconstruct::dump_coefficients(string file)
{
  ofstream out(file.c_str());
  if (!out) {
    cerr << "Reconstruct::dump_coefficients: can't open <"
	 << file << "> . . . " << endl;
    return;
  }
  out.setf(ios::scientific);
  out.precision(8);

  double fac;

  // Header
  out << "# Cosine coefficients only" << endl;
  out << "#" << setw(4) << "l" << setw(5) << "m";
  for (int i=1; i<=nmax; i++) {
    ostringstream ostr;
    ostr << "n=" << i;
    out << setw(18) << ostr.str();
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


void profile_header(ostream& out)
{
  const int ncol = 6, siz = 15;
  const char *lab[] = {"Position", "Dens(x)", "Potl(x)", 
		       "Dens(y)", "Potl(y)", "Dpot(x)"};
  const char *sep[] = {"#", "+", "+", "+", "+", "+"};
  const char *cnt[] = {"[1]", "[2]", "[3]", "[4]", "[5]", "[6]"};

  out << left << setfill('-');
  for (int i=0; i<ncol; i++) out << sep[i] << setw(siz-1) << '-';
  out << endl << setfill(' ');
  for (int i=0; i<ncol; i++) out << setw(2) << sep[i] << setw(siz-2) << lab[i];
  out << endl;
  for (int i=0; i<ncol; i++) out << setw(2) << sep[i] << setw(siz-2) << cnt[i];
  out << endl << setfill('-');
  for (int i=0; i<ncol; i++) out << sep[i] << setw(siz-1) << '-';
  out << endl << setfill(' ');
}


int 
main(int argc, char** argv)
{
  //====================
  // Parse command line
  //====================

  try {
    if (config.parse_args(argc, argv)) return -1;
    assign_values();
  }
  catch (const char *msg) {
    cerr << msg << endl;
    return -1;
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

    profile_header(cout);
    profile_header(out);

    double x, dx = 2.0*xmax/(numx-1);
    double xm, xp, dp;
    for (int i=0; i<numx; i++) {
      x = -xmax + dx*i;
      xm = max<double>(rmin, x*(1.0-delta));
      xp = min<double>(rmax, x*(1.0+delta));
      dp = (
	    recon.potential_eval(xp, 0.0, 0.0) -
	    recon.potential_eval(xm, 0.0, 0.0)
	    ) / (xp - xm);
      cout << setw(15) << x
	   << setw(15) << recon.density_eval(x, 0.0, 0.0)
	   << setw(15) << recon.potential_eval(x, 0.0, 0.0)
	   << setw(15) << recon.density_eval(delr, x, 0.0)
	   << setw(15) << recon.potential_eval(delr, x, 0.0)
	   << setw(15) << dp
	   << endl;
      out  << setw(15) << x
	   << setw(15) << recon.density_eval(x, 0.0, 0.0)
	   << setw(15) << recon.potential_eval(x, 0.0, 0.0)
	   << setw(15) << recon.density_eval(delr, x, 0.0)
	   << setw(15) << recon.potential_eval(delr, x, 0.0)
	   << setw(15) << dp
	   << endl;
    }

    if (surface) {
      double x, y, dxy = 2.0*xmax/(numx-1);
      float z;

      ofstream out(string(outfile + ".surf").c_str());

      if (out) {
	for (int j=0; j<numx; j++) {
	  y = -xmax + dxy*j;
	  for (int i=0; i<numx; i++) {
	    x = -xmax + dxy*i;
	  
	    out << setw(15) << x << setw(15) << y
		<< setw(15) << recon.potential_eval(x, y, 0.0)
		<< setw(15) << recon.density_eval(x, y, 0.0) << endl;
	  }
	  out << endl;
	}

      } else {
	cerr << "Could not open <" << outfile + ".surf" << ">" << endl;
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
