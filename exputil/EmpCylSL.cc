#include <filesystem>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <memory>
#include <limits>
#include <string>

#include <yaml-cpp/yaml.h>

#include <Progress.H>		// Progress bar
#include <TableGrid.H>		// For three-dimensional array

#include <interp.H>
#include <Timer.H>
#include <thread>
#include <exp_thread.h>
#include <EXPException.H>

#include <Eigen/Eigenvalues>

#include <config.h>
#ifdef HAVE_VTK
#include <VtkPCA.H>
#endif

#ifdef HAVE_OMP_H
#include <omp.h>		// For multithreading basis construction
#endif

#include <numerical.H>
#include <gaussQ.H>
#include <EmpCylSL.H>
#include <DataGrid.H>

#include <libvars.H>
using namespace __EXP__;	// For reference to n-body globals

#undef  TINY
#define TINY 1.0e-16


bool     EmpCylSL::DENS            = false;
bool     EmpCylSL::PCAVAR          = false;
bool     EmpCylSL::PCAVTK          = false;
bool     EmpCylSL::PCAEOF          = false;
bool     EmpCylSL::PCADRY          = true;
bool     EmpCylSL::logarithmic     = false;
bool     EmpCylSL::enforce_limits  = false;
int      EmpCylSL::CMAPR           = 1;
int      EmpCylSL::CMAPZ           = 1;
int      EmpCylSL::NUMX            = 256;
int      EmpCylSL::NUMY            = 128;
int      EmpCylSL::NOUT            = 12;
int      EmpCylSL::NUMR            = 2000;
unsigned EmpCylSL::VFLAG           = 0;
unsigned EmpCylSL::VTKFRQ          = 1;
double   EmpCylSL::HEXP            = 1.0;
double   EmpCylSL::RMIN            = 0.001;
double   EmpCylSL::RMAX            = 20.0;
double   EmpCylSL::HFAC            = 0.2;
double   EmpCylSL::PPOW            = 4.0;
bool     EmpCylSL::NewCache        = true;
bool     EmpCylSL::NewCoefs        = true;
 

EmpCylSL::EmpModel EmpCylSL::mtype = Exponential;

std::map<EmpCylSL::EmpModel, std::string> EmpCylSL::EmpModelLabs =
  { {Exponential, "Exponential"},
    {Gaussian,    "Gaussian"   },
    {Plummer,     "Plummer"    },
    {Power,       "Power"      },
    {Deproject,   "Deproject"  }
  };


EmpCylSL::EmpCylSL()
{
  NORDER     = 0;
  coefs_made = std::vector<short>(multistep+1, false);
  eof_made   = false;
  defSampT   = 1;
  sampT      = 1;
  tk_type    = None;
  EVEN_M     = false;
  MMIN       = 0;
  MLIM       = std::numeric_limits<int>::max();
  NMIN       = 0;
  NLIM       = std::numeric_limits<int>::max();
  EvenOdd    = false;
  Neven      = 0;
  Nodd       = 0;
  minSNR     = std::numeric_limits<double>::max();
  maxSNR     = 0.0;
  cachefile  = default_cache;

  if (DENS)
    MPItable = 4;
  else
    MPItable = 3;

  cylmass = 0.0;
  cylmass_made = false;
  cylmass1 = std::vector<double>(nthrds);

  hallfile = "";
}

EmpCylSL::~EmpCylSL(void)
{
  // Nothing . . . 
}

EmpCylSL::EmpCylSL(int nmax, int lmax, int mmax, int nord, 
		   double ascale, double hscale, int nodd,
		   std::string cachename)
{
  // Use default name?
  if (cachename.size()) cachefile = cachename;
  else                  cachefile = default_cache;

  // Sanity check
  if (lmax <= mmax) {
    if (myid==0) {
      std::cout << "EmpCylSL: lmax must be greater than mmax for consistency"
		<< std::endl
		<< "EmpCylSL: setting lmax=" << mmax + 1
		<< " but you probably want lmax >> mmax"
		<< std::endl;
    }
    lmax = mmax + 1;
  }

  NMAX      = nmax;
  MMAX      = mmax;
  LMAX      = lmax;
  NORDER    = nord;
  MMIN      = 0;
  MLIM      = std::numeric_limits<int>::max();
  NMIN      = 0;
  NLIM      = std::numeric_limits<int>::max();
  EvenOdd   = false;
  Neven     = 0;
  Nodd      = 0;


  ASCALE   = ascale;
  HSCALE   = hscale;
  pfac     = 1.0/sqrt(ascale);
  ffac     = pfac/ascale;
  dfac     = ffac/ascale;

  EVEN_M   = false;

  // Set number of even and odd terms
  //
  if (nodd>=0 and nodd<=NORDER) {
    EvenOdd  = true;
    Neven    = NORDER - nodd;
    Nodd     = nodd;
  }

  // Enable MPI code for more than one node
  //
  if (numprocs>1) SLGridSph::mpi = 1;

  if (DENS)
    MPItable = 4;
  else
    MPItable = 3;

  coefs_made = std::vector<short>(multistep+1, false);
  eof_made   = false;

  sampT        = 1;
  defSampT     = 1;
  tk_type      = None;

  cylmass      = 0.0;
  cylmass1     = std::vector<double>(nthrds);
  cylmass_made = false;

  hallfile     = "";
  minSNR       = std::numeric_limits<double>::max();
  maxSNR       = 0.0;
}


void EmpCylSL::reset(int numr, int lmax, int mmax, int nord, 
		     double ascale, double hscale, int nodd,
		     std::string cachename)
{
  // Reset cache file name
  if (cachename.size()) cachefile = cachename;

  // Option sanity check
  if (lmax <= mmax) {
    if (myid==0) {
      std::cout << "EmpCylSL: lmax must be greater than mmax for consistency"
		<< std::endl
		<< "EmpCylSL: setting lmax=" << mmax + 1
		<< " but you probably want lmax >> mmax"
		<< std::endl;
    }
    lmax = mmax + 1;
  }

  NMAX     = numr;
  MMAX     = mmax;
  LMAX     = lmax;
  NORDER   = nord;
  MMIN     = 0;
  MLIM     = std::numeric_limits<int>::max();
  NMIN     = 0;
  NLIM     = std::numeric_limits<int>::max();
  EvenOdd  = false;

  // Set number of even and odd terms
  //
  if (nodd>=0 and nodd<=NORDER) {
    EvenOdd  = true;
    Neven    = NORDER - nodd;
    Nodd     = nodd;
  }

  ASCALE = ascale;
  HSCALE = hscale;
  pfac   = 1.0/sqrt(ascale);
  ffac   = pfac/ascale;
  dfac   = ffac/ascale;

  SLGridSph::mpi = 1;		// Turn on MPI
  ortho = std::make_shared<SLGridSph>(make_sl(), LMAX, NMAX, NUMR,
				      RMIN, RMAX*0.99, false, 1, 1.0);

  coefs_made = std::vector<short>(multistep+1, false);
  eof_made = false;

  if (DENS)
    MPItable = 4;
  else
    MPItable = 3;

  sampT = 1;
  defSampT = 1;

  cylmass = 0.0;
  cylmass1.resize(nthrds);
  cylmass_made = false;

  minSNR  = std::numeric_limits<double>::max();
  maxSNR  = 0.0;
}

void EmpCylSL::cacheInfo(const std::string& cachefile)
{
  auto node = getHeader(cachefile);

  std::cout << std::string(60, '-') << std::endl
	    << "Cache parameters for EmpCylSL: " << cachefile << std::endl
	    << std::string(60, '-') << std::endl;
  for (YAML::const_iterator it=node.begin(); it!=node.end(); ++it) {
    std::cout << std::left << std::setw(20) << it->first.as<std::string>()
	      << ": " << it->second.as<std::string>() << std::endl;
  }
  std::cout << std::string(60, '-') << std::endl;
}

void EmpCylSL::create_deprojection(double H, double Rf, int NUMR, int NINT,
				   AxiDiskPtr func)
{
  LegeQuad lq(NINT);

  std::vector<double> rr(NUMR), rl(NUMR), sigI(NUMR), rhoI(NUMR, 0.0);

  double Rmin = log(RMIN);
  double Rmax = log(RMAX);

  double dr = (Rmax - Rmin)/(NUMR-1);

  // Compute surface mass density, Sigma(R)
  //
  for (int i=0; i<NUMR; i++) {
    double r = Rmin + dr*i;

    // Save for finite difference
    //
    rl[i] = r;
    r = exp(r);
    rr[i] = r;

    // Interval by Legendre
    //
    sigI[i] = 0.0;
    for (int n=0; n<NINT; n++) {
      double y   = lq.knot(n);
      double y12 = 1.0 - y*y;
      double z   = y/sqrt(y12)*H;

      sigI[i] += lq.weight(n)*2.0*H*pow(y12, -1.5)*(*func)(r*Rf, z);
    }
  }

  Linear1d surf(rl, sigI);
  
  // Now, compute Abel inverion integral
  //
  for (int i=0; i<NUMR; i++) {
    double r = rr[i];

    // Interval by Legendre
    //
    rhoI[i] = 0.0;
    for (int n=0; n<NINT; n++) {
      double x   = lq.knot(n);
      double x12 = 1.0 - x*x;
      double z   = x/sqrt(x12);
      double R   = sqrt(z*z + r*r);
      double lR  = log(R);

      rhoI[i]   += lq.weight(n)*2.0*pow(x12, -1.5)*surf.eval(lR);
    }
  }

  std::vector<double> rho(NUMR), mass(NUMR);

  Linear1d intgr(rl, rhoI);

  for (int i=0; i<NUMR; i++)
    rho[i] = -intgr.deriv(rl[i])/(2.0*M_PI*rr[i]*rr[i]);

  mass[0] = 0.0;
  for (int i=1; i<NUMR; i++) {
    double rlst = rr[i-1], rcur = rr[i];
    mass[i] = mass[i-1] + 2.0*M_PI*(rlst*rlst*rho[i-1] + rcur*rcur*rho[i])*(rcur - rlst);
  }

  // Debug
  //
  if (VFLAG & 1) {
    std::ostringstream outf; outf << "deproject_sl." << myid;
    std::ofstream out(outf.str());
    if (out) {
      for (int i=0; i<NUMR; i++)
	out << std::setw(18) << rl[i]
	    << std::setw(18) << rr[i]
	    << std::setw(18) << rho[i]
	    << std::setw(18) << mass[i]
	    << std::endl;
    }
  }

  // Finalize
  //
  densRg = Linear1d(rl, rho);
  massRg = Linear1d(rl, mass);
  mtype  = Deproject;
}


/*
  Note that the produced by the following three routines
  are in dimensionless units
*/
double EmpCylSL::massR(double R)
{
  double ans=0.0, fac, arg;

  switch (mtype) {
  case Exponential:
    ans = 1.0 - (1.0 + R)*exp(-R); 
    break;
  case Gaussian:
    arg = 0.5*R*R;
    ans = 1.0 - exp(-arg);
    break;
  case Plummer:
    fac = R/(1.0+R);
    ans = pow(fac, 3.0);
    break;
  case Power:
    {
      double z  = R + 1.0;
      double a1 = PPOW - 1.0;
      double a2 = PPOW - 2.0;
      double a3 = PPOW - 3.0;
      ans =  0.5*a1*a2*a3 * ( (1.0 - pow(z, -a3))/a3 -
			      (1.0 - pow(z, -a2))/a2 * 2.0 +
			      (1.0 - pow(z, -a1))/a1 );
    }
    break;
  case Deproject:
    if (R < RMIN) ans = 0.0;
    else if (R>=RMAX) ans = massRg.eval(RMAX);
    else ans = massRg.eval(log(R));
    break;
  }

  return ans;
}

double EmpCylSL::densR(double R)
{
  double ans=0.0, fac, arg;

  switch (mtype) {
  case Exponential:
    ans = exp(-R)/(4.0*M_PI*R);
    break;
  case Gaussian:
    arg = 0.5*R*R;
    ans = exp(-arg)/(4.0*M_PI*R);
    break;
  case Plummer:
    fac = 1.0/(1.0+R);
    ans = 3.0*pow(fac, 4.0)/(4.0*M_PI);
    break;
  case Power:
    {
      double z  = R + 1.0;
      double a1 = PPOW - 1.0;
      double a2 = PPOW - 2.0;
      double a3 = PPOW - 3.0;
      ans =  0.125*a1*a2*a3/M_PI * pow(z, -PPOW);
    }
    break;
  case Deproject:
    if (R < RMIN) ans = densRg.eval(RMIN);
    else if (R>=RMAX) ans = 0.0;
    else ans = densRg.eval(log(R));
    break;
  }

  return ans;
}

SphModTblPtr EmpCylSL::make_sl()
{
  const int number = 10000;

  r.resize(number);
  d.resize(number);
  m.resize(number);
  p.resize(number);

  std::vector<double> mm(number);
  std::vector<double> pw(number);

				// ------------------------------------------
				// Debug sanity check
				// ------------------------------------------
  if (myid==0) {
    std::cout << "EmpCylSL::make_sl(): making SLGridSph with <"
	      << EmpModelLabs[mtype] << "> model" << std::endl;
  }

				// ------------------------------------------
				// Make radial, density and mass array
				// ------------------------------------------
  double dr;
  if (logarithmic)
    dr = (log(RMAX) - log(RMIN))/(number - 1);
  else
    dr = (RMAX - RMIN)/(number - 1);

  for (int i=0; i<number; i++) {
    if (logarithmic)
      r[i] = RMIN*exp(dr*i);
    else
      r[i] = RMIN + dr*i;

    m[i] = massR(r[i]);
    d[i] = densR(r[i]);
  }

  mm[0] = 0.0;
  pw[0] = 0.0;
  for (int i=1; i<number; i++) {
    mm[i] = mm[i-1] +
      2.0*M_PI*(r[i-1]*r[i-1]*d[i-1] + r[i]*r[i]*d[i])*
      (r[i] - r[i-1]);
    pw[i] = pw[i-1] +
      2.0*M_PI*(r[i-1]*d[i-1] + r[i]*d[i])*(r[i] - r[i-1]);
  }

  for (int i=0; i<number; i++) 
    p[i] = -mm[i]/(r[i]+1.0e-10) - (pw[number-1] - pw[i]);

  if (VFLAG & 16) {
    ostringstream outf;
    outf << "test_adddisk_sl." << myid;
    ofstream out(outf.str().c_str());
    for (int i=0; i<number; i++) {
      out 
	<< setw(15) << r[i] 
	<< setw(15) << d[i] 
	<< setw(15) << m[i] 
	<< setw(15) << p[i] 
	<< setw(15) << mm[i] 
	<< endl;
    }
    out.close();
  }

  return std::make_shared<SphericalModelTable>
    (number, r.data(), d.data(), m.data(), p.data());
}


void EmpCylSL::send_eof_grid()
{

#ifdef SEND_DEBUG
  auto blab = [](auto& s, auto& t, auto id, auto m, auto v) {
		 std::cout << "[" << id << "] " << s << "  " << t
			     << " (" << m << ", " << v << ")" << std::endl;
		 };
  
  std::cout << "[" << myid << "] (MMAX, rank3)=("
	    << MMAX << ", " << rank3 << ")" << std::endl;
#else
  auto blab = [](auto& s, auto& t, auto id, auto m, auto v) {};
#endif

  // Send to slaves
  // 
  for (int m=0; m<=MMAX; m++) {

    // Grids in X--Y
    // 
    for (int v=0; v<rank3; v++) {

      blab("before", "potC", myid, m, v);
      MPI_Bcast(potC[m][v].data(), potC[m][v].size(),
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
      blab("after", "potC", myid, m, v);

      blab("before", "rforceC", myid, m, v);
      MPI_Bcast(rforceC[m][v].data(), rforceC[m][v].size(),
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
      blab("after", "rforceC", myid, m, v);

      blab("before", "zforceC", myid, m, v);
      MPI_Bcast(zforceC[m][v].data(), zforceC[m][v].size(),
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
      blab("after", "zforceC", myid, m, v);

      if (DENS) {
	blab("before", "densC", myid, m, v);
	MPI_Bcast(densC[m][v].data(), densC[m][v].size(),
		  MPI_DOUBLE, 0, MPI_COMM_WORLD);
	blab("after", "densC", myid, m, v);
      }
    }
  }

  for (int m=1; m<=MMAX; m++) {
    
    // Grids in X--Y
    //
    for (int v=0; v<rank3; v++) {

      blab("before", "potS", myid, m, v);
      MPI_Bcast(potS[m][v].data(), potS[m][v].size(),
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
      blab("after", "potS", myid, m, v);

      blab("before", "rforceS", myid, m, v);
      MPI_Bcast(rforceS[m][v].data(), rforceS[m][v].size(),
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
      blab("after", "rforceS", myid, m, v);

      blab("before", "zforceS", myid, m, v);
      MPI_Bcast(zforceS[m][v].data(), zforceS[m][v].size(),
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
      blab("after", "zforceS", myid, m, v);

      if (DENS) {
	blab("before", "densS", myid, m, v);
	MPI_Bcast(densS[m][v].data(), densS[m][v].size(),
		    MPI_DOUBLE, 0, MPI_COMM_WORLD);
	blab("after", "densS", myid, m, v);
      }
      
    }
      
  }
}


int EmpCylSL::read_eof_header(const std::string& eof_file)
{
  std::ifstream in(eof_file.c_str());
  if (!in) {
    std::cerr << "---- EmpCylSL::read_eof_header: error opening file named <" 
	      << eof_file << ">" << std::endl;
    return 0;
  }

  // Attempt to read magic number
  //
  unsigned int tmagic;
  in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

  if (tmagic == hmagic) {

      // YAML size
      //
      unsigned ssize;
      in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));

      // Make and read char buffer
      //
      auto buf = std::make_unique<char[]>(ssize+1);
      in.read(buf.get(), ssize);
      buf[ssize] = 0;		// Null terminate

      YAML::Node node;
      
      try {
	node = YAML::Load(buf.get());
      }
      catch (YAML::Exception& error) {
	if (myid==0)
	  std::cerr << "YAML: error parsing <" << buf.get() << "> "
		    << "in " << __FILE__ << ":" << __LINE__ << std::endl
		    << "YAML error: " << error.what() << std::endl;
	throw error;
      }

      // Get parameters
      //
      MMAX   = node["mmax"  ].as<int>();
      NUMX   = node["numx"  ].as<int>();
      NUMY   = node["numy"  ].as<int>();
      NMAX   = node["nmax"  ].as<int>();
      NORDER = node["norder"].as<int>();
      DENS   = node["dens"  ].as<bool>();
      RMIN   = node["rmin"  ].as<double>();
      RMAX   = node["rmax"  ].as<double>();
      ASCALE = node["ascl"  ].as<double>();
      HSCALE = node["hscl"  ].as<double>();

      if (node["cmap"])		// Backwards compatibility
	CMAPR  = node["cmap"  ].as<int>();
      else
	CMAPR  = node["cmapr" ].as<int>();

      if (node["cmapz"])	// Backwards compatibility
	CMAPZ  = node["cmapz" ].as<int>();


  } else {

    // Rewind file
    //
    in.clear();
    in.seekg(0);

    int tmp;

    in.read((char *)&MMAX,   sizeof(int));
    in.read((char *)&NUMX,   sizeof(int));
    in.read((char *)&NUMY,   sizeof(int));
    in.read((char *)&NMAX,   sizeof(int));
    in.read((char *)&NORDER, sizeof(int));
    in.read((char *)&tmp,    sizeof(int)); 
    if (tmp) DENS = true; else DENS = false;
    in.read((char *)&CMAPR,  sizeof(int)); 
    in.read((char *)&RMIN,   sizeof(double));
    in.read((char *)&RMAX,   sizeof(double));
    in.read((char *)&ASCALE, sizeof(double));
    in.read((char *)&HSCALE, sizeof(double));
  }
  
  if (myid==0) {
    cout << setfill('-') << setw(70) << '-' << endl;
    cout << " Cylindrical parameters read from <" << eof_file << ">" << endl;
    cout << setw(70) << '-' << endl;
    cout << "MMAX="   << MMAX   << endl;
    cout << "NUMX="   << NUMX   << endl;
    cout << "NUMY="   << NUMY   << endl;
    cout << "NMAX="   << NMAX   << endl;
    cout << "NORDER=" << NORDER << endl;
    cout << "DENS="   << DENS   << endl;
    cout << "CMAPR="  << CMAPR  << endl;
    cout << "CMAPZ="  << CMAPZ  << endl;
    cout << "RMIN="   << RMIN   << endl;
    cout << "RMAX="   << RMAX   << endl;
    cout << "ASCALE=" << ASCALE << endl;
    cout << "HSCALE=" << HSCALE << endl;
    cout << setw(70) << '-' << endl << setfill(' ');
  }

  return 1;
}

int EmpCylSL::read_eof_file(const string& eof_file)
{
  read_eof_header(eof_file);

  pfac = 1.0/sqrt(ASCALE);
  ffac = pfac/ASCALE;
  dfac = ffac/ASCALE;

  SLGridSph::mpi = 1;		// Turn on MPI
  ortho = std::make_shared<SLGridSph>(make_sl(), LMAX, NMAX, NUMR,
					RMIN, RMAX*0.99, false, 1, 1.0);

  setup_eof();
  setup_accumulation();

				// Master tries to read table
  int retcode;
  if (myid==0) retcode = cache_grid(0, eof_file);
  MPI_Bcast(&retcode, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (!retcode) return 0;
				// Send table to slave processes
  send_eof_grid();

  if (myid==0) 
    std::cout << "---- EmpCylSL::read_cache: table forwarded to all processes"
	      << std::endl;


  eof_made = true;
  coefs_made = std::vector<short>(multistep+1, false);

  return 1;
}

int EmpCylSL::read_cache(void)
{
  setup_table();
  setup_accumulation();

				// Master tries to read table
  int retcode;
  if (myid==0) retcode = cache_grid(0);
  MPI_Bcast(&retcode, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (!retcode) return 0;
				// Send table to slave processes
  send_eof_grid();

  if (myid==0) 
    std::cout << "---- EmpCylSL::read_cache: table forwarded to all processes"
	      << std::endl;


  eof_made = true;
  coefs_made = std::vector<short>(multistep+1, false);

  return 1;
}


template <typename U>
std::string compare_out(std::string str, U one, U two)
{
  std::ostringstream sout;

  std::transform(str.begin(), str.end(),str.begin(), ::toupper);

  sout << std::setw(15) << std::right << str
       << std::setw(15) << std::right << one
       << std::setw(15) << std::right << two
       << std::endl;

  return sout.str();
}

int EmpCylSL::cache_grid(int readwrite, string cachename)
{

  // Option to reset cache file name
  if (cachename.size()) cachefile = cachename;

  if (readwrite) {
    
    if (std::filesystem::exists(cachefile)) {
      std::cerr << "---- EmpCylSL::cache_grid: cache file <"
		<< cachefile << "> exists" << std::endl;
      try {
	std::filesystem::rename(cachefile, cachefile + ".bak");
      }
      catch(std::filesystem::filesystem_error const& ex) {
        std::cout << "---- EmpCylSL::cache_grid write error: "
		  << "what():  " << ex.what()  << std::endl
		  << "path1(): " << ex.path1() << std::endl
		  << "path2(): " << ex.path2() << std::endl;
	MPI_Abort(MPI_COMM_WORLD, 12);
	return 0;
      }

      std::cout << "---- EmpCylSL::cache_grid: existing file backed up to <"
		<< cachefile + ".bak>" << std::endl;
    }

    std::ofstream out(cachefile);
    if (!out) {
      std::cerr << "---- EmpCylSL::cache_grid: error opening file for writing"
		<< std::endl;
      return 0;
    }

    if (NewCache) {

      // This is a node of simple {key: value} pairs.  More general
      // content can be added as needed.
      YAML::Node node;

      node["model" ] = EmpModelLabs[mtype];
      node["mmax"  ] = MMAX;
      node["numx"  ] = NUMX;
      node["numy"  ] = NUMY;
      node["nmax"  ] = NMAX;
      node["norder"] = NORDER;
      node["neven" ] = Neven;
      node["nodd"  ] = Nodd;
      node["dens"  ] = DENS;
      node["cmapr" ] = CMAPR;
      node["cmapz" ] = CMAPZ;
      node["rmin"  ] = RMIN;
      node["rmax"  ] = RMAX;
      node["ascl"  ] = ASCALE;
      node["hscl"  ] = HSCALE;
      node["cmass" ] = cylmass;

      // Serialize the node
      //
      YAML::Emitter y; y << node;

      // Get the size of the string
      //
      unsigned int hsize = strlen(y.c_str());

      // Write magic #
      //
      out.write(reinterpret_cast<const char *>(&hmagic),   sizeof(unsigned int));

      // Write YAML string size
      //
      out.write(reinterpret_cast<const char *>(&hsize),    sizeof(unsigned int));

      // Write YAML string
      //
      out.write(reinterpret_cast<const char *>(y.c_str()), hsize);

    } else {

      // Old-style header
      //
      const int one  = 1;
      const int zero = 0;
      double time    = 0.0;

      out.write((const char *)&MMAX,    sizeof(int));
      out.write((const char *)&NUMX,    sizeof(int));
      out.write((const char *)&NUMY,    sizeof(int));
      out.write((const char *)&NMAX,    sizeof(int));
      out.write((const char *)&NORDER,  sizeof(int));
      if (DENS) out.write((const char *)&one,  sizeof(int));
      else      out.write((const char *)&zero, sizeof(int));
      out.write((const char *)&CMAPR,   sizeof(int));
      out.write((const char *)&RMIN,    sizeof(double));
      out.write((const char *)&RMAX,    sizeof(double));
      out.write((const char *)&ASCALE,  sizeof(double));
      out.write((const char *)&HSCALE,  sizeof(double));
      out.write((const char *)&cylmass, sizeof(double));
      out.write((const char *)&time,    sizeof(double));
    }

    // Write table
    //
    for (int m=0; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&potC[m][v](ix, iy), sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&rforceC[m][v](ix, iy), sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&zforceC[m][v](ix, iy), sizeof(double));
	  
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      out.write((const char *)&densC[m][v](ix, iy), sizeof(double));

	}
      }
    }

    for (int m=1; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&potS[m][v](ix, iy), sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&rforceS[m][v](ix, iy), sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    out.write((const char *)&zforceS[m][v](ix, iy), sizeof(double));
	
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      out.write((const char *)&densS[m][v](ix, iy), sizeof(double));
	}
      }

    }
    
  }
  else {

    std::ifstream in(cachefile.c_str());
    if (!in) {
      cerr << "---- EmpCylSL::cache_grid: error opening file" << endl;
      return 0;
    }

    int mmax, numx, numy, nmax, norder, tmp, cmapr, cmapz;
    bool dens=false;
    double rmin, rmax, ascl, hscl, time;


    // Attempt to read magic number
    //
    unsigned int tmagic;
    in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

    if (tmagic == hmagic) {
      // YAML size
      //
      unsigned ssize;
      in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));

      // Make and read char buffer
      //
      auto buf = std::make_unique<char[]>(ssize+1);
      in.read(buf.get(), ssize);
      buf[ssize] = 0;		// Null terminate

      YAML::Node node;
      
      try {
	node = YAML::Load(buf.get());
      }
      catch (YAML::Exception& error) {
	if (myid)
	  std::cerr << "YAML: error parsing <" << buf.get() << "> "
		    << "in " << __FILE__ << ":" << __LINE__ << std::endl
		    << "YAML error: " << error.what() << std::endl;
	throw error;
      }

      // Get parameters
      //
      mmax    = node["mmax"  ].as<int>();
      numx    = node["numx"  ].as<int>();
      numy    = node["numy"  ].as<int>();
      nmax    = node["nmax"  ].as<int>();
      norder  = node["norder"].as<int>();
      dens    = node["dens"  ].as<bool>();
      rmin    = node["rmin"  ].as<double>();
      rmax    = node["rmax"  ].as<double>();
      ascl    = node["ascl"  ].as<double>();
      hscl    = node["hscl"  ].as<double>();
      cylmass = node["cmass" ].as<double>();
      
      if (node["time"]) 	// Backwards compatibility; 'time' is optional
	time    = node["time"  ].as<double>();
      else
	time    = 0.0;

      if (node["cmap"]) 	// Backwards compatibility; parameter change
	cmapr   = node["cmap" ].as<int>();
      else
	cmapr   = node["cmapr"].as<int>();

      if (node["cmapz"])	// Backwards compatibility; parameter change
	cmapz = node["cmapz"].as<int>();
      else
	cmapz = CMAPZ;

    } else {
				// Rewind file
      in.clear();
      in.seekg(0);

      in.read((char *)&mmax,    sizeof(int));
      in.read((char *)&numx,    sizeof(int));
      in.read((char *)&numy,    sizeof(int));
      in.read((char *)&nmax,    sizeof(int));
      in.read((char *)&norder,  sizeof(int));
      in.read((char *)&tmp,     sizeof(int));    if (tmp) dens = true;
      in.read((char *)&cmapr,   sizeof(int));
      in.read((char *)&rmin,    sizeof(double));
      in.read((char *)&rmax,    sizeof(double));
      in.read((char *)&ascl,    sizeof(double));
      in.read((char *)&hscl,    sizeof(double));
      in.read((char *)&cylmass, sizeof(double));
      in.read((char *)&time,    sizeof(double));

      cmapz = 1;
    }

				// Spot check compatibility
    if ( (MMAX    != mmax   ) |
	 (NUMX    != numx   ) |
	 (NUMY    != numy   ) |
	 (NMAX    != nmax   ) |
	 (NORDER  != norder ) |
	 (DENS    != dens   ) |
	 (CMAPR   != cmapr  ) |
	 (fabs(rmin-RMIN)>1.0e-12 ) |
	 (fabs(rmax-RMAX)>1.0e-12 ) |
	 (fabs(ascl-ASCALE)>1.0e-12 ) |
	 (fabs(hscl-HSCALE)>1.0e-12 )
	 ) 
      {
	cout << std::setw(15) << std::right << "Key"
	     << std::setw(15) << std::right << "Wanted"
	     << std::setw(15) << std::right << "Cached"
	     << std::endl
	     << std::setw(15) << std::right << "------"
	     << std::setw(15) << std::right << "------"
	     << std::setw(15) << std::right << "------"
	     << std::endl;

	cout << compare_out("mmax",   MMAX,   mmax);
	cout << compare_out("numx",   NUMX,   numx);
	cout << compare_out("numy",   NUMY,   numy);
	cout << compare_out("nmax",   NMAX,   nmax);
	cout << compare_out("norder", NORDER, norder);
	cout << compare_out("dens",   DENS,   dens);
	cout << compare_out("cmapr",  CMAPR,  cmapr);
	cout << compare_out("cmapz",  CMAPZ,  cmapz);
	cout << compare_out("rmin",   RMIN,   rmin);
	cout << compare_out("rmax",   RMAX,   rmax);
	cout << compare_out("ascale", ASCALE, ascl);
	cout << compare_out("hscale", HSCALE, hscl);

	return 0;
      }

				// Read table

    Eigen::MatrixXd work(NUMY+1, NUMX+1);

    for (int m=0; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&potC[m][v](ix, iy), sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&rforceC[m][v](ix, iy), sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&zforceC[m][v](ix, iy), sizeof(double));
	  
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      in.read((char *)&densC[m][v](ix, iy), sizeof(double));

	}
      }

    }

    for (int m=1; m<=MMAX; m++) {

      for (int v=0; v<rank3; v++) {

	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&potS[m][v](ix, iy), sizeof(double));
	
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&rforceS[m][v](ix, iy), sizeof(double));
	  
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++)
	    in.read((char *)&zforceS[m][v](ix, iy), sizeof(double));
	
	if (DENS) {
	  for (int ix=0; ix<=NUMX; ix++)
	    for (int iy=0; iy<=NUMY; iy++)
	      in.read((char *)&densS[m][v](ix, iy), sizeof(double));
	}
      }

    }
    
    Rtable = M_SQRT1_2 * RMAX;
    XMIN = r_to_xi(RMIN*ASCALE);
    XMAX = r_to_xi(Rtable*ASCALE);
    dX = (XMAX - XMIN)/NUMX;
    
    YMIN = z_to_y(-Rtable*ASCALE);
    YMAX = z_to_y( Rtable*ASCALE);
    dY = (YMAX - YMIN)/NUMY;
    
    std::cout << "---- EmpCylSL::cache_grid: file read successfully"
	      << std::endl;
  }

  return 1;
}

YAML::Node EmpCylSL::getHeader(const std::string& cachefile)
{
  YAML::Node node;

  std::ifstream in(cachefile);
  if (!in) {
    std::ostringstream sout;
    sout << "EmpCylSL::getHeader: could not open cache file <" << cachefile << ">";
    std::runtime_error(sout.str());
  }

  // Attempt to read magic number
  //
  unsigned int tmagic;
  in.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

  if (tmagic == hmagic) {

    // YAML size
    //
    unsigned ssize;
    in.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));

    // Make and read char buffer
    //
    auto buf = std::make_unique<char[]>(ssize+1);
    in.read(buf.get(), ssize);
    buf[ssize] = 0;		// Null terminate
    
    try {
      node = YAML::Load(buf.get());
    }
    catch (YAML::Exception& error) {
      if (myid==0)
	std::cerr << "YAML: error parsing <" << buf.get() << "> "
		  << "in " << __FILE__ << ":" << __LINE__ << std::endl
		  << "YAML error: " << error.what() << std::endl;
      throw error;
    }
  } else {
    std::ostringstream sout;
    sout << "EmpCylSL::getHeader: invalid cache file <" << cachefile << ">";
    std::runtime_error(sout.str());
  }

  return node;
}
    

void EmpCylSL::receive_eof(int request_id, int MM)
{
  int type, icnt, off;
  int mm;

  MPI_Recv(&type, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, 
	   MPI_COMM_WORLD, &status);

  int current_source = status.MPI_SOURCE;

  if (VFLAG & 16)
    cerr << "Master beginning to receive from " << current_source 
	 << " . . . " << endl;

  MPI_Recv(&mm, 1, MPI_INT, current_source, MPI_ANY_TAG, 
	   MPI_COMM_WORLD, &status);

  if (VFLAG & 16)
    cerr << "Master receiving from " << current_source << ": type=" << type 
	 << "   M=" << mm << endl;

				// Receive rest of data

  for (int n=0; n<NORDER; n++) {
    MPI_Recv(&mpi_double_buf2[MPIbufsz*(MPItable*n+0)], 
	     MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+1, 
	     MPI_COMM_WORLD, &status);

    MPI_Recv(&mpi_double_buf2[MPIbufsz*(MPItable*n+1)], 
	     MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+2, 
	     MPI_COMM_WORLD, &status);

    MPI_Recv(&mpi_double_buf2[MPIbufsz*(MPItable*n+2)], 
	     MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+3, 
	     MPI_COMM_WORLD, &status);
    if (DENS)
      MPI_Recv(&mpi_double_buf2[MPIbufsz*(MPItable*n+3)], 
	       MPIbufsz, MPI_DOUBLE, current_source, 13 + MPItable*n+4, 
	       MPI_COMM_WORLD, &status);
  }
  

				// Send slave new orders
  if (request_id >=0) {
    MPI_Send(&request_id, 1, MPI_INT, current_source, 1, MPI_COMM_WORLD);
    MPI_Send(&MM, 1, MPI_INT, current_source, 2, MPI_COMM_WORLD);
  }
  else {
    MPI_Send(&request_id, 1, MPI_INT, current_source, 1, MPI_COMM_WORLD);
  }
  
				// Read from buffers

  for (int n=0; n<NORDER; n++) {
  
    off = MPIbufsz*(MPItable*n+0);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	potC[mm][n](ix, iy)  = mpi_double_buf2[off+icnt++];
      else
	potS[mm][n](ix, iy)  = mpi_double_buf2[off+icnt++];
	


    off = MPIbufsz*(MPItable*n+1);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	rforceC[mm][n](ix, iy)  = mpi_double_buf2[off+icnt++];
      else
	rforceS[mm][n](ix, iy)  = mpi_double_buf2[off+icnt++];
	

    off = MPIbufsz*(MPItable*n+2);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
      if (type)
	zforceC[mm][n](ix, iy)  = mpi_double_buf2[off+icnt++];
      else
	zforceS[mm][n](ix, iy)  = mpi_double_buf2[off+icnt++];
    

    if (DENS) {
      off = MPIbufsz*(MPItable*n+3);
      icnt = 0;
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++)
	  if (type)
	    densC[mm][n](ix, iy)  = mpi_double_buf2[off+icnt++];
	  else
	    densS[mm][n](ix, iy)  = mpi_double_buf2[off+icnt++];
    }
  }
  
  if (VFLAG & 16)
    cerr << "Master finished receiving: type=" << type << "   M=" 
	 << mm << endl;

  return;

}

void EmpCylSL::compute_eof_grid(int request_id, int m)
{
  // Check for existence of ortho and create if necessary
  //
  if (not ortho)
    ortho = std::make_shared<SLGridSph>(make_sl(), LMAX, NMAX, NUMR,
					  RMIN, RMAX*0.99, false, 1, 1.0);


  //  Read in coefficient matrix or
  //  make grid if needed
  double fac1, fac2, dens, potl, potr, pott, fac3, fac4;
  
  int icnt, off;
  
  for (int v=0; v<NORDER; v++) {
    tpot[v].setZero();
    trforce[v].setZero();
    tzforce[v].setZero();
    if (DENS) tdens[v].setZero();
  }

  for (int ix=0; ix<=NUMX; ix++) {

    double x = XMIN + dX*ix;
    double r = xi_to_r(x);

    for (int iy=0; iy<=NUMY; iy++) {

      double y = YMIN + dY*iy;
      double z = y_to_z(y);

      double rr = sqrt(r*r + z*z) + 1.0e-18;

      ortho->get_pot(potd, rr/ASCALE);
      ortho->get_force(dpot, rr/ASCALE);
      if (DENS) ortho->get_dens(dend, rr/ASCALE);

      double costh = z/rr;
      dlegendre_R(LMAX, costh, legs[0], dlegs[0]);
      
      for (int v=0; v<NORDER; v++) {

	for (int ir=0; ir<NMAX; ir++) {

	  for (int l=m; l<=LMAX; l++) {

	    fac1 = 1.;//sqrt((2.0*l+1.0)/(4.0*M_PI));

	    if (m==0) {
	      fac2 = legs[0](l, m);

	      dens = fac2*dend(l, ir) * dfac;
	      potl = fac2*potd(l, ir) * pfac;
	      potr = fac2*dpot(l, ir) * ffac;
	      pott = dlegs[0](l, m)*potd(l, ir) * pfac;

	    } else {

	      fac2 = M_SQRT2;
	      fac3 = fac2 * legs[0](l, m);
	      fac4 = fac2 * dlegs[0](l, m);
	      
	      dens = fac3*dend(l, ir) * dfac;
	      potl = fac3*potd(l, ir) * pfac;
	      potr = fac3*dpot(l, ir) * ffac;
	      pott = fac4*potd(l, ir) * pfac;
	    }
	    
	    int nn = ir + NMAX*(l-m);

	    tpot[v](ix, iy) +=  ef(nn, v) * potl;

	    trforce[v](ix, iy) += 
	      -ef(nn, v) * (potr*r/rr - pott*z*r/(rr*rr*rr));

	    tzforce[v](ix, iy) += 
	      -ef(nn, v) * (potr*z/rr + pott*r*r/(rr*rr*rr));

	    if (DENS) 
	      tdens[v](ix, iy) +=  ef(nn, v) * dens;
	  }
	}
      }
    }
  }

				// Send stuff back to master
      
  MPI_Send(&request_id, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
  MPI_Send(&m, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
      
    
  for (int n=0; n<NORDER; n++) {

				// normalization factors
    if (DENS)
      tdens[n] *= 0.25/M_PI;

				// Potential

    off = MPIbufsz*(MPItable*n+0);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = tpot[n](ix, iy);
    
    if (VFLAG & 16)
      cerr << "Slave " << setw(4) << myid << ": with request_id=" << request_id
	   << ", M=" << m << " send Potential" << endl;

    MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	     13 + MPItable*n+1, MPI_COMM_WORLD);

				// R force

    off = MPIbufsz*(MPItable*n+1);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = trforce[n](ix, iy);
    
    if (VFLAG & 16)
      cerr << "Slave " << setw(4) << myid << ": with request_id=" << request_id
	   << ", M=" << m << " sending R force" << endl;

    MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	     13 + MPItable*n+2, MPI_COMM_WORLD);

				// Z force

    off = MPIbufsz*(MPItable*n+2);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = tzforce[n](ix, iy);
    
    if (VFLAG & 16)
      cerr << "Slave " << setw(4) << myid << ": with request_id=" << request_id
	   << ", M=" << m << " sending Z force" << endl;


    MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	     13 + MPItable*n+3, MPI_COMM_WORLD);

				// Density

    if (DENS) {
      off = MPIbufsz*(MPItable*n+3);
      icnt = 0;
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++)
	  mpi_double_buf2[off + icnt++] = tdens[n](ix, iy);
    
      if (VFLAG & 16)
	cerr << "Slave " << setw(4) << myid 
	     << ": with request_id=" << request_id
	     << ", M=" << m << " sending Density" << endl;

      MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	       13 + MPItable*n+4, MPI_COMM_WORLD);

    }

  }

}

void EmpCylSL::compute_even_odd(int request_id, int m)
{
  // check for ortho
  //
  if (not ortho)
    ortho = std::make_shared<SLGridSph>(make_sl(),
					  LMAX, NMAX, NUMR, RMIN, RMAX*0.99,
					  false, 1, 1.0);

  double fac1, fac2, dens, potl, potr, pott, fac3, fac4;
  
  int icnt, off;
  
  for (int v=0; v<NORDER; v++) {
    tpot[v].setZero();
    trforce[v].setZero();
    tzforce[v].setZero();
    if (DENS) tdens[v].setZero();
  }

  for (int ix=0; ix<=NUMX; ix++) {

    double x = XMIN + dX*ix;
    double r = xi_to_r(x);

    for (int iy=0; iy<=NUMY; iy++) {

      double y = YMIN + dY*iy;
      double z = y_to_z(y);

      double rr = sqrt(r*r + z*z) + 1.0e-18;

      ortho->get_pot(potd, rr/ASCALE);
      ortho->get_force(dpot, rr/ASCALE);
      if (DENS) ortho->get_dens(dend, rr/ASCALE);

      double costh = z/rr;
      dlegendre_R(LMAX, costh, legs[0], dlegs[0]);
      
      // Do the even eigenfunction first
      //
      for (int v=0; v<Neven; v++) {

	for (int ir=0; ir<NMAX; ir++) {

	  for (int il=0; il<lE[m].size(); il++) {

	    int l = lE[m][il];	// Even l values first
	    
	    fac1 = 1.;//sqrt((2.0*l+1.0)/(4.0*M_PI));

	    if (m==0) {
	      fac2 = legs[0](l, m);

	      dens = fac2*dend(l, ir) * dfac;
	      potl = fac2*potd(l, ir) * pfac;
	      potr = fac2*dpot(l, ir) * ffac;
	      pott = dlegs[0](l, m)*potd(l, ir) * pfac;

	    } else {

	      fac2 = M_SQRT2;
	      fac3 = fac2 * legs[0](l, m);
	      fac4 = fac2 * dlegs[0](l, m);
	      
	      dens = fac3*dend(l, ir) * dfac;
	      potl = fac3*potd(l, ir) * pfac;
	      potr = fac3*dpot(l, ir) * ffac;
	      pott = fac4*potd(l, ir) * pfac;
	    }
	    
	    int nn = ir + NMAX*il;

	    tpot[v](ix, iy) +=  efE(nn, v) * potl;

	    trforce[v](ix, iy) += 
	      -efE(nn, v) * (potr*r/rr - pott*z*r/(rr*rr*rr));

	    tzforce[v](ix, iy) += 
	      -efE(nn, v) * (potr*z/rr + pott*r*r/(rr*rr*rr));

	    if (DENS) 
	      tdens[v](ix, iy) +=  efE(nn, v) * dens;
	  }
	}
      }

      for (int v=Neven; v<NORDER; v++) {

	int w = v - Neven;	// Index in odd eigenfunctions

	for (int ir=0; ir<NMAX; ir++) {

	  for (int il=0; il<lO[m].size(); il++) {

	    int l = lO[m][il];

	    fac1 = 1.;//sqrt((2.0*l+1.0)/(4.0*M_PI));

	    if (m==0) {
	      fac2 = legs[0](l, m);

	      dens = fac2*dend(l, ir) * dfac;
	      potl = fac2*potd(l, ir) * pfac;
	      potr = fac2*dpot(l, ir) * ffac;
	      pott = dlegs[0](l, m)*potd(l, ir) * pfac;

	    } else {

	      fac2 = M_SQRT2;
	      fac3 = fac2 * legs[0](l, m);
	      fac4 = fac2 * dlegs[0](l, m);
	      
	      dens = fac3*dend(l, ir) * dfac;
	      potl = fac3*potd(l, ir) * pfac;
	      potr = fac3*dpot(l, ir) * ffac;
	      pott = fac4*potd(l, ir) * pfac;
	    }
	    
	    int nn = ir + NMAX*il;

	    tpot[v](ix, iy) +=  efO(nn, w) * potl;

	    trforce[v](ix, iy) += 
	      -efO(nn, w) * (potr*r/rr - pott*z*r/(rr*rr*rr));

	    tzforce[v](ix, iy) += 
	      -efO(nn, w) * (potr*z/rr + pott*r*r/(rr*rr*rr));

	    if (DENS) 
	      tdens[v](ix, iy) +=  efO(nn, w) * dens;
	  }
	}
      }
    }
  }
  
				// Send stuff back to master
      
  MPI_Send(&request_id, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
  MPI_Send(&m, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
      
    
  for (int n=0; n<NORDER; n++) {

				// normalization factors
    if (DENS)
      tdens[n] *= 0.25/M_PI;

				// Potential

    off = MPIbufsz*(MPItable*n+0);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = tpot[n](ix, iy);
    
    if (VFLAG & 16)
      cerr << "Slave " << setw(4) << myid << ": with request_id=" << request_id
	   << ", M=" << m << " send Potential" << endl;

    MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	     13 + MPItable*n+1, MPI_COMM_WORLD);

				// R force

    off = MPIbufsz*(MPItable*n+1);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = trforce[n](ix, iy);
    
    if (VFLAG & 16)
      cerr << "Slave " << setw(4) << myid << ": with request_id=" << request_id
	   << ", M=" << m << " sending R force" << endl;

    MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	     13 + MPItable*n+2, MPI_COMM_WORLD);

				// Z force

    off = MPIbufsz*(MPItable*n+2);
    icnt = 0;
    for (int ix=0; ix<=NUMX; ix++)
      for (int iy=0; iy<=NUMY; iy++)
	mpi_double_buf2[off + icnt++] = tzforce[n](ix, iy);
    
    if (VFLAG & 16)
      cerr << "Slave " << setw(4) << myid << ": with request_id=" << request_id
	   << ", M=" << m << " sending Z force" << endl;


    MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	     13 + MPItable*n+3, MPI_COMM_WORLD);

				// Density

    if (DENS) {
      off = MPIbufsz*(MPItable*n+3);
      icnt = 0;
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++)
	  mpi_double_buf2[off + icnt++] = tdens[n](ix, iy);
    
      if (VFLAG & 16)
	cerr << "Slave " << setw(4) << myid 
	     << ": with request_id=" << request_id
	     << ", M=" << m << " sending Density" << endl;

      MPI_Send(&mpi_double_buf2[off], MPIbufsz, MPI_DOUBLE, 0, 
	       13 + MPItable*n+4, MPI_COMM_WORLD);

    }

  }

}


void EmpCylSL::setup_accumulation(int mlevel)
{
  if (accum_cos.size()==0) {	// First time only

    accum_cos.resize(MMAX+1);
    accum_sin.resize(MMAX+1);

    cosL.resize(multistep+1);
    cosN.resize(multistep+1);
    sinL.resize(multistep+1);
    sinN.resize(multistep+1);

    howmany1.resize(multistep+1);
    howmany .resize(multistep+1, 0);

    for (unsigned M=0; M<=multistep; M++) {
      cosL[M] = std::make_shared<VectorD2>(nthrds);
      cosN[M] = std::make_shared<VectorD2>(nthrds);
      sinL[M] = std::make_shared<VectorD2>(nthrds);
      sinN[M] = std::make_shared<VectorD2>(nthrds);
      
      for (int nth=0; nth<nthrds; nth++) {
	cosL(M)[nth].resize(MMAX+1);
	cosN(M)[nth].resize(MMAX+1);
	sinL(M)[nth].resize(MMAX+1);
	sinN(M)[nth].resize(MMAX+1);
      }

      howmany1[M].resize(nthrds, 0);
    }

    if (VFLAG & 16)
      cerr << "Slave " << setw(4) << myid 
	   << ": tables allocated, MMAX=" << MMAX << endl;

    cylmass_made = false;

    for (unsigned M=0; M<=multistep; M++) {
      
      for (int nth=0; nth<nthrds; nth++) {
	
	for (int m=0; m<=MMAX; m++) {
	  
	  cosN(M)[nth][m].resize(NORDER);
	  cosL(M)[nth][m].resize(NORDER);
	  
	  if (m>0) {
	    sinN(M)[nth][m].resize(NORDER);
	    sinL(M)[nth][m].resize(NORDER);
	  }
	}
      }
    }
    
    for (int m=0; m<=MMAX; m++) {
      accum_cos[m].resize(NORDER);
      if (m>0) accum_sin[m].resize(NORDER);
    }
    
    for (unsigned M=0; M<=multistep; M++) {
      for (int nth=0; nth<nthrds; nth++) {
	for (int m=0; m<=MMAX; m++) {
	  cosN(M)[nth][m].setZero();
	  if (m>0) sinN(M)[nth][m].setZero();
	}
      }
    }

    if (PCAVAR and sampT>0) {
      for (int nth=0; nth<nthrds; nth++) {
	for (unsigned T=0; T<sampT; T++) {
	  numbT1[nth][T] = 0;
	  massT1[nth][T] = 0.0;
	  covV[nth][T].resize(MMAX+1);
	  covM[nth][T].resize(MMAX+1);
	  for (int mm=0; mm<=MMAX; mm++) {
	    covV[nth][T][mm].resize(NORDER);
	    covM[nth][T][mm].resize(NORDER, NORDER);
	  }
	}
      }
    }
  }

  // Zero values on every pass
  //
  for (int m=0; m<=MMAX; m++) {
    accum_cos[m].setZero();
    if (m>0) accum_sin[m].setZero();
  }

  if ( (PCAVAR or PCAEOF) and mlevel==0 and sampT>0) {

    for (int nth=0; nth<nthrds; nth++) {

      if (PCAEOF) {
	for (auto & v : tvar[nth]) v.setZero();
      }

      if (PCAVAR) {
	minSNR = std::numeric_limits<double>::max();
	maxSNR = 0.0;

	for (unsigned T=0; T<sampT; T++) {
	  numbT1[nth][T] = 0;
	  massT1[nth][T] = 0.0;
	  for (int mm=0; mm<=MMAX; mm++) {
	    covV[nth][T][mm].setZero();
	    covM[nth][T][mm].setZero();
	  }
	}
      }
    }
  }

  // Reset particle counter
  //
  howmany[mlevel] = 0;

  // Swap buffers
  //
  auto p  = cosL[mlevel];
  cosL[mlevel] = cosN[mlevel];
  cosN[mlevel] = p;
    
  p       = sinL[mlevel];
  sinL[mlevel] = sinN[mlevel];
  sinN[mlevel] = p;
    
  // Clean current coefficient files
  //
  for (int nth=0; nth<nthrds; nth++) {
    
    howmany1[mlevel][nth] = 0;
      
    for (int m=0; m<=MMAX; m++) {
      cosN(mlevel)[nth][m].setZero();
      if (m>0) sinN(mlevel)[nth][m].setZero();
    }
  }
    
  coefs_made[mlevel] = false;

  // DONE
}

void EmpCylSL::init_pca()
{
  if (PCAVAR or PCAEOF) {
    if (PCAVAR) {
      if (defSampT) sampT = defSampT;
      else          sampT = floor(sqrt(nbodstot));
    }

    pthread_mutex_init(&used_lock, NULL);

    if (PCAEOF)
      tvar.resize(nthrds);

    if (PCAVAR) {
      covV  .resize(nthrds);
      covM  .resize(nthrds);
      numbT1.resize(nthrds);
      massT1.resize(nthrds);
      numbT .resize(sampT, 0);
      massT .resize(sampT, 0);
    }

    for (int nth=0; nth<nthrds;nth++) {
      if (PCAEOF) {
	tvar[nth].resize(MMAX + 1);
	for (auto & v : tvar[nth]) v.resize(rank3, rank3);
      }

      if (PCAVAR) {
	minSNR = std::numeric_limits<double>::max();
	maxSNR = 0.0;

	numbT1[nth].resize(sampT, 0);
	massT1[nth].resize(sampT, 0);

	covV[nth].resize(sampT);
	covM[nth].resize(sampT);
	for (unsigned T=0; T<sampT; T++) {
	  covV[nth][T].resize(MMAX+1);
	  covM[nth][T].resize(MMAX+1);
	  for (int mm=0; mm<=MMAX; mm++) {
	    covV[nth][T][mm].resize(rank3);
	    covM[nth][T][mm].resize(rank3, rank3);
	  }
	}
      }
    }
  }
}

void EmpCylSL::setup_table()
{
  // Create storage for EOF tables
  //
  rank2   = NMAX*(LMAX+1);
  rank3   = NORDER;
    
  Rtable  = M_SQRT1_2 * RMAX;
  XMIN    = r_to_xi(RMIN*ASCALE);
  XMAX    = r_to_xi(Rtable*ASCALE);
  dX      = (XMAX - XMIN)/NUMX;

  YMIN    = z_to_y(-Rtable*ASCALE);
  YMAX    = z_to_y( Rtable*ASCALE);
  dY      = (YMAX - YMIN)/NUMY;

  potC    .resize(MMAX+1);
  rforceC .resize(MMAX+1);
  zforceC .resize(MMAX+1);
  
  potS    .resize(MMAX+1);
  rforceS .resize(MMAX+1);
  zforceS .resize(MMAX+1);
  
  if (DENS) {
    densC .resize(MMAX+1);
    densS .resize(MMAX+1);
  }

  for (int m=0; m<=MMAX; m++) {

    potC[m]   .resize(rank3);
    rforceC[m].resize(rank3);
    zforceC[m].resize(rank3);
    if (DENS) densC[m].resize(rank3);
    
    for (int v=0; v<rank3; v++) {
      potC   [m][v].resize(NUMX+1, NUMY+1);
      rforceC[m][v].resize(NUMX+1, NUMY+1);
      zforceC[m][v].resize(NUMX+1, NUMY+1);
      if (DENS) densC[m][v].resize(NUMX+1, NUMY+1);
    }
  }
  
  for (int m=1; m<=MMAX; m++) {
    
    potS[m]   .resize(rank3);
    rforceS[m].resize(rank3);
    zforceS[m].resize(rank3);
    if (DENS) densS[m].resize(rank3);
    
    for (int v=0; v<rank3; v++) {
      potS   [m][v].resize(NUMX+1, NUMY+1);
      rforceS[m][v].resize(NUMX+1, NUMY+1);
      zforceS[m][v].resize(NUMX+1, NUMY+1);
      if (DENS) densS[m][v].resize(NUMX+1, NUMY+1);
    }
  }

  vc.resize(nthrds);
  vs.resize(nthrds);
  for (int i=0; i<nthrds; i++) {
    vc[i].resize(max<int>(1,MMAX)+1, rank3);
    vs[i].resize(max<int>(1,MMAX)+1, rank3);
  }

}

void EmpCylSL::setup_eof()
{
  if (SC.size()==0 and SCe.size()==0) {

    setup_table();

    tpot   .resize(NORDER);
    trforce.resize(NORDER);
    tzforce.resize(NORDER);
    if (DENS) tdens.resize(NORDER);

    for (int n=0; n<NORDER; n++) {
      tpot[n].resize(NUMX+1, NUMY+1);
      trforce[n].resize(NUMX+1, NUMY+1);
      tzforce[n].resize(NUMX+1, NUMY+1);
      if (DENS) tdens[n].resize(NUMX+1, NUMY+1);
    }

    if (EvenOdd) {
      SCe.resize(nthrds);
      SSe.resize(nthrds);

      SCo.resize(nthrds);
      SSo.resize(nthrds);

      for (int nth=0; nth<nthrds; nth++) {
	SCe[nth].resize(MMAX+1);
	SSe[nth].resize(MMAX+1);

	SCo[nth].resize(MMAX+1);
	SSo[nth].resize(MMAX+1);
      }

      lE.resize(MMAX+1);
      lO.resize(MMAX+1);

    } else {
      SC.resize(nthrds);
      SS.resize(nthrds);

      for (int nth=0; nth<nthrds; nth++) {
	SC[nth].resize(MMAX+1);
	SS[nth].resize(MMAX+1);
      }
    }

    potd.resize(LMAX+1, NMAX);
    dpot.resize(LMAX+1, NMAX);
    dend.resize(LMAX+1, NMAX);

    cosm .resize(nthrds);
    sinm .resize(nthrds);
    legs .resize(nthrds);
    dlegs.resize(nthrds);
    for (int i=0; i<nthrds; i++) {
      cosm[i].resize(LMAX+1);
      sinm[i].resize(LMAX+1);
      legs[i].resize(LMAX+1, LMAX+1);
      dlegs[i].resize(LMAX+1, LMAX+1);
    }

    if (EvenOdd) {
      varE.resize(MMAX+1);
      varO.resize(MMAX+1);

      for (int m=0; m<=MMAX; m++) {

	lE[m].clear();	// Store even l values
	lO[m].clear();	// Store odd  l values

	for (int l=m; l<=LMAX; l++) {
	  // Symmetry for the associated Legendre polynomials
	  //   |   depends on whether l+m is even or odd
	  //   v        
	  if ((l+m) % 2==0) lE[m].push_back(l);
	  else              lO[m].push_back(l);
	}
	
	int Esiz = lE[m].size();
	int Osiz = lO[m].size();

	varE[m].resize(NMAX*Esiz, NMAX*Esiz);
	varO[m].resize(NMAX*Osiz, NMAX*Osiz);
      }

    } else {
      var.resize(MMAX+1);
      for (int m=0; m<=MMAX; m++)
	var[m].resize(NMAX*(LMAX-m+1), NMAX*(LMAX-m+1));
    }

    for (int nth=0; nth<nthrds; nth++) {

      for (int m=0; m<=MMAX; m++) {

	if (EvenOdd) {
	  int Esiz = lE[m].size();
	  int Osiz = lO[m].size();

	  SCe[nth][m].resize(NMAX*Esiz);
	  SCo[nth][m].resize(NMAX*Osiz);

	  if (m) {
	    SSe[nth][m].resize(NMAX*Esiz);
	    SSo[nth][m].resize(NMAX*Osiz);
	  }

	  for (int i=0; i<NMAX*Esiz; i++) {
	    SCe[nth][m][i].resize(NMAX*Esiz);
	    if (m) SSe[nth][m][i].resize(NMAX*Esiz);
	  }

	  for (int i=0; i<NMAX*Osiz; i++) {
	    SCo[nth][m][i].resize(NMAX*Osiz);
	    if (m) SSo[nth][m][i].resize(NMAX*Osiz);
	  }

	} else {

	  SC[nth][m].resize(NMAX*(LMAX-m+1));
	  if (m) SS[nth][m].resize(NMAX*(LMAX-m+1));

	  for (int i=0; i<NMAX*(LMAX-m+1); i++) {

	    SC[nth][m][i].resize(NMAX*(LMAX-m+1));
	    if (m) SS[nth][m][i].resize(NMAX*(LMAX-m+1));
	  
	  }
	}
      }
    
    }

    table.resize(nthrds);
    facC .resize(nthrds);
    facS .resize(nthrds);
    for (int i=0; i<nthrds; i++) {
      table[i].resize(LMAX+1, NMAX);
      facC [i].resize(NMAX, LMAX+1);
      facS [i].resize(NMAX, LMAX+1);
    }

    MPIbufsz = (NUMX+1)*(NUMY+1);

    mpi_double_buf2.resize(MPIbufsz*NORDER*MPItable);
    mpi_double_buf3.resize(rank3);
  }

  for (int nth=0; nth<nthrds; nth++) {
    for (int m=0; m<=MMAX; m++)  {
      
      if (EvenOdd) {
	int Esiz = lE[m].size();
	int Osiz = lO[m].size();

	for (int i=0; i<NMAX*Esiz; i++)  {
	  for (int j=0; j<NMAX*Esiz; j++)  {
	    SCe[nth][m][i][j] = 0.0;
	    if (m>0) SSe[nth][m][i][j] = 0.0;
	  }
	}

	for (int i=0; i<NMAX*Osiz; i++)  {
	  for (int j=0; j<NMAX*Osiz; j++)  {
	    SCo[nth][m][i][j] = 0.0;
	    if (m>0) SSo[nth][m][i][j] = 0.0;
	  }
	}

      } else {
	for (int i=0; i<NMAX*(LMAX-m+1); i++)  {
	  for (int j=0; j<NMAX*(LMAX-m+1); j++)  {
	    SC[nth][m][i][j] = 0.0;
	    if (m>0) SS[nth][m][i][j] = 0.0;
	  }
	}
      }
    }
  }

  eof_made = false;
}


// Create EOF from target density and spherical basis
//
void EmpCylSL::generate_eof(int numr, int nump, int numt, 
			    double (*func)
			    (double R, double z, double phi, int M) )
{
  Timer timer;
  if (VFLAG & 16) timer.start();

  // Create spherical orthogonal basis if necessary
  //
  if (not ortho)
    ortho = std::make_shared<SLGridSph>(make_sl(), LMAX, NMAX, NUMR,
					  RMIN, RMAX*0.99, false, 1, 1.0);
  // Initialize fixed variables and storage
  //
  setup_eof();

  LegeQuad lr(numr);
  LegeQuad lt(numt);
  double dphi = 2.0*M_PI/nump;


#ifdef HAVE_OMP_H
  omp_set_dynamic(0);		// Explicitly disable dynamic teams
  omp_set_num_threads(nthrds);	// OpenMP set up
#endif

  std::shared_ptr<progress::progress_display> progress;
  if (VFLAG & 16 && myid==0) {
    std::cout << std::endl << "Quadrature loop progress" << std::endl;
    progress = std::make_shared<progress::progress_display>(numr);
  }

  int cntr = 0;			// Loop counter for spreading load to nodes
  
  // *** Radial quadrature loop
  //
  for (int qr=0; qr<numr; qr++) { 

    // Diagnostic timing output for MPI process loop
    //
    if (VFLAG & 16 && myid==0) {
      ++(*progress);
    }    

    if (cntr++ % numprocs != myid) continue;

    double xi = XMIN + (XMAX - XMIN) * lr.knot(qr);
    double rr = xi_to_r(xi);
    ortho->get_pot(table[0], rr/ASCALE);

    // *** cos(theta) quadrature loop
    //
#pragma omp parallel for
    for (int qt=0; qt<numt; qt++) {
#ifdef HAVE_OMP_H
      int id = omp_get_thread_num();
#else
      int id = 0;
#endif
	
      double costh = -1.0 + 2.0*lt.knot(qt);
      double R     = rr * sqrt(1.0 - costh*costh);
      double z     = rr * costh;
      
      legendre_R(LMAX, costh, legs[id]);

      double jfac = dphi*2.0*lt.weight(qt)*(XMAX - XMIN)*lr.weight(qr) 
	* rr*rr / d_xi_to_r(xi);
      
      // *** Phi quadrature loop
      //
      for (int qp=0; qp<nump; qp++) {

	double phi = dphi*qp;
	sinecosine_R(LMAX, phi, cosm[id], sinm[id]);

	// *** m loop
	//
	for (int m=0; m<=MMAX; m++) {

	  // Get the target density for this position and azimuthal index
	  //
	  double dens = (*func)(R, z, phi, m) * jfac;

	  // *** ir loop
	  //
	  for (int ir=0; ir<NMAX; ir++) {

	    // *** l loop
	    //
	    for (int l=m; l<=LMAX; l++) {
		
	      double ylm = pfac * legs[id](l, m);

	      if (m==0) {

		facC[id](ir, l-m) = ylm*table[0](l, ir);
		 
                // very deep debug for model problems. Not necessarily worth a debug flag.
		//if (std::isnan(facC[id](ir, l-m))) cerr << "Process " << myid << ": EmpCylSL. Has nan in facC. table=[" << table[0](l, ir) << "] ylm=[" << ylm << "]" << endl;
	      }
	      else {
		
		if (nump==1) {
		  facC[id](ir, l-m) = ylm*table[0](l, ir)*0.5;
		  facS[id](ir, l-m) = ylm*table[0](l, ir)*0.5;
		} else {
		  facC[id](ir, l-m) = ylm*table[0](l, ir)*cosm[id][m];
		  facS[id](ir, l-m) = ylm*table[0](l, ir)*sinm[id][m];
		}
	      }
	      
	    } // *** l loop
	    
	  } // *** ir loop

	  for (int ir1=0; ir1<NMAX; ir1++) {
	    
	    if (EvenOdd) {
	      
	      // Even l loop
	      //
	      for (int il1=0; il1<lE[m].size(); il1++) {

		int l1  = lE[m][il1];
		int nn1 = ir1 + NMAX*il1;

		if (m==0) {
		
		  for (int ir2=0; ir2<NMAX; ir2++) {

		    for (int il2=0; il2<lE[m].size(); il2++) {

		      int l2  = lE[m][il2];
		      int nn2 = ir2 + NMAX*il2;
		    
		      SCe[id][m][nn1][nn2] += 
			facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * dens;
		    }
		  }
		
		} else {
		
		  for (int ir2=0; ir2<NMAX; ir2++) {

		    for (int il2=0; il2<lE[m].size(); il2++) {
		      
		      int l2  = lE[m][il2];
		      int nn2 = ir2 + NMAX*il2;

		      SCe[id][m][nn1][nn2] += 
			facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * dens;
		      
		      SSe[id][m][nn1][nn2] += 
			facS[id](ir1, l1-m)*facS[id](ir2, l2-m) * dens;
		    }
		  }
		}
	      
	      } // *** Even l loop
	    
	      
	      // Odd l loop
	      //
	      for (int il1=0; il1<lO[m].size(); il1++) {

		int l1  = lO[m][il1];
		int nn1 = ir1 + NMAX*il1;

		if (m==0) {
		
		  for (int ir2=0; ir2<NMAX; ir2++) {

		    for (int il2=0; il2<lO[m].size(); il2++) {

		      int l2  = lO[m][il2];
		      int nn2 = ir2 + NMAX*il2;
		    
		      SCo[id][m][nn1][nn2] += 
			facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * dens;
		    }
		  }
		
		} else {
		
		  for (int ir2=0; ir2<NMAX; ir2++) {

		    for (int il2=0; il2<lO[m].size(); il2++) {

		      int l2  = lO[m][il2];
		      int nn2 = ir2 + NMAX*il2;
		      
		      SCo[id][m][nn1][nn2] += 
			facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * dens;
		    
		      SSo[id][m][nn1][nn2] += 
			facS[id](ir1, l1-m)*facS[id](ir2, l2-m) * dens;
		    }
		  }
		}
	      
	      } // *** Odd l loop
	      

	    } // *** END: EvenOdd==true block
	    else {

	      for (int l1=m; l1<=LMAX; l1++) {

		int nn1 = ir1 + NMAX*(l1-m);

		if (m==0) {
		
		  for (int ir2=0; ir2<NMAX; ir2++) {
		    for (int l2=m; l2<=LMAX; l2++) {
		      int nn2 = ir2 + NMAX*(l2-m);
		    
		      SC[id][m][nn1][nn2] += 
			facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * dens;
		    }
		  }
		
		} else {
		
		  for (int ir2=0; ir2<NMAX; ir2++) {

		    for (int l2=m; l2<=LMAX; l2++) {
		      
		      int nn2 = ir2 + NMAX*(l2-m);
		      
		      SC[id][m][nn1][nn2] += 
			facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * dens;
		      
		      SS[id][m][nn1][nn2] += 
			facS[id](ir1, l1-m)*facS[id](ir2, l2-m) * dens;
		    }
		  }
		}
		
	      } // *** l loop
	    
	    } // *** EvenOdd block

	  } // *** ir loop
	  
	} // *** m loop

      } // *** phi quadrature loop

    } // *** cos(theta) quadrature loop

  } // *** r quadrature loop
  
  if (VFLAG & 16) {
    auto t = timer.stop();
    if (myid==0) {
      std::cout << std::endl
		<< std::setw( 8) << "Process"
		<< std::setw(18) << "Completion time"
		<< std::endl
		<< std::setw( 8) << "-------"
		<< std::setw(18) << "---------------"
		<< std::endl;
    }

    for (int n=0; n<numprocs; n++) {
      if (n==myid) {
	std::cout << std::setw( 8) << myid
		  << std::setw(18) << t
		  << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    if (myid==0) {
      std::cout << std::endl;
    }

    timer.reset();
    timer.start();
  }

  //
  // Now, we are ready to make the EOF basis
  //

  make_eof();

  if (VFLAG & 16) {
    cout << "Process " << setw(4) << myid << ": completed basis in " 
	 << timer.stop() << " seconds"
	 << endl;
  }

  //
  // We still need to make the coefficients
  //

  coefs_made = std::vector<short>(multistep+1, false);

}


void EmpCylSL::accumulate_eof(double r, double z, double phi, double mass, 
			      int id, int mlevel)
{
  if (not ortho)
    ortho = std::make_shared<SLGridSph>
      (make_sl(), LMAX, NMAX, NUMR, RMIN, RMAX*0.99, false, 1, 1.0);
  if (eof_made) {
    if (VFLAG & 2)
      cerr << "accumulate_eof: Process " << setw(4) << myid << ", Thread " 
	   << id << " calling setup_eof()" << endl;
    setup_eof();
  }

  double rr = sqrt(r*r + z*z);

  if (rr/ASCALE>Rtable) return;

  ortho->get_pot(table[id], rr/ASCALE);
  double costh = z/(rr+1.0e-18);
  legendre_R(LMAX, costh, legs[id]);
  sinecosine_R(LMAX, phi, cosm[id], sinm[id]);

  int nn1, nn2;

  // *** m loop
  for (int m=0; m<=MMAX; m++) {

    // *** ir loop
    for (int ir=0; ir<NMAX; ir++) {

      // *** l loop
      for (int l=m; l<=LMAX; l++) {

	double ylm = pfac * legs[id](l, m);

	if (m==0) {

	  facC[id](ir, l-m) = ylm*table[id](l, ir);

	}
	else {

	  facC[id](ir, l-m) = ylm*table[id](l, ir)*cosm[id][m];
	  facS[id](ir, l-m) = ylm*table[id](l, ir)*sinm[id][m];

	}

      } // *** l loop

    } // *** ir loop

    for (int ir1=0; ir1<NMAX; ir1++) {

      if (EvenOdd) {
	int Esiz = lE[m].size();
	int Osiz = lO[m].size();

	for (int il1=0; il1<Esiz; il1++) {

	  int l1 = lE[m][il1];
	  nn1 = ir1 + NMAX*il1;

	  if (m==0) {
	  
	    for (int ir2=0; ir2<NMAX; ir2++) {

	      for (int il2=0; il2<Esiz; il2++) {

		int l2 = lE[m][il2];
		nn2 = ir2 + NMAX*il2;

		SCe[id][m][nn1][nn2] += 
		  facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * mass;
	      }
	    }

	  } else {

	    for (int ir2=0; ir2<NMAX; ir2++) {

	      for (int il2=0; il2<Esiz; il2++) {

		int l2 = lE[m][il2];
		nn2 = ir2 + NMAX*il2;
		
		SCe[id][m][nn1][nn2] += 
		  facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * mass;
		
		SSe[id][m][nn1][nn2] += 
		  facS[id](ir1, l1-m)*facS[id](ir2, l2-m) * mass;
	      }
	    }
	  }

	} // *** Even l loop

	for (int il1=0; il1<Osiz; il1++) {

	  int l1 = lO[m][il1];
	  nn1 = ir1 + NMAX*il1;

	  if (m==0) {
	  
	    for (int ir2=0; ir2<NMAX; ir2++) {

	      for (int il2=0; il2<Osiz; il2++) {

		int l2 = lO[m][il2];
		nn2 = ir2 + NMAX*il2;

		SCo[id][m][nn1][nn2] += 
		  facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * mass;
	      }
	    }

	  } else {

	    for (int ir2=0; ir2<NMAX; ir2++) {

	      for (int il2=0; il2<Osiz; il2++) {

		int l2 = lO[m][il2];
		nn2 = ir2 + NMAX*il2;

		SCo[id][m][nn1][nn2] += 
		  facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * mass;
		
		SSo[id][m][nn1][nn2] += 
		  facS[id](ir1, l1-m)*facS[id](ir2, l2-m) * mass;
	      }
	    }
	  }

	} // *** Odd l loop

      } else {
      
	for (int l1=m; l1<=LMAX; l1++) {
	  nn1 = ir1 + NMAX*(l1-m);

	  if (m==0) {
	  
	    for (int ir2=0; ir2<NMAX; ir2++) {
	      for (int l2=m; l2<=LMAX; l2++) {
		nn2 = ir2 + NMAX*(l2-m);
		
		SC[id][m][nn1][nn2] += 
		  facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * mass;
	      }
	    }
	    
	  } else {

	    for (int ir2=0; ir2<NMAX; ir2++) {
	      for (int l2=m; l2<=LMAX; l2++) {
		nn2 = ir2 + NMAX*(l2-m);
		
		SC[id][m][nn1][nn2] += 
		  facC[id](ir1, l1-m)*facC[id](ir2, l2-m) * mass;

		SS[id][m][nn1][nn2] += 
		  facS[id](ir1, l1-m)*facS[id](ir2, l2-m) * mass;
	      }
	    }
	  }
	  
	} // *** l loop
	
      } // *** EvenOdd loop

    } // *** ir loop

  } // *** m loop
  
}

void EmpCylSL::make_eof(void)
{
  Timer timer;
  int icnt;
  double tmp;

  if (MPIin_eof.size()==0) {
    MPIin_eof .resize(rank2*(rank2+1)/2);
    MPIout_eof.resize(rank2*(rank2+1)/2);
  }
  
  //
  //  Sum up over threads
  //

  for (int nth=1; nth<nthrds; nth++) {

    for (int mm=0; mm<=MMAX; mm++) {

      if (EvenOdd) {
	int Esiz = lE[mm].size();
	int Osiz = lO[mm].size();

	for (int i=0; i<NMAX*Esiz; i++)
	  for (int j=i; j<NMAX*Esiz; j++)
	    SCe[0][mm][i][j] += SCe[nth][mm][i][j];

	for (int i=0; i<NMAX*Osiz; i++)
	  for (int j=i; j<NMAX*Osiz; j++)
	    SCo[0][mm][i][j] += SCo[nth][mm][i][j];

      } else {
	for (int i=0; i<NMAX*(LMAX-mm+1); i++)
	  for (int j=i; j<NMAX*(LMAX-mm+1); j++)
	    SC[0][mm][i][j] += SC[nth][mm][i][j];
      }

    }

    for (int mm=1; mm<=MMAX; mm++) {

      if (EvenOdd) {
	int Esiz = lE[mm].size();
	int Osiz = lO[mm].size();

	for (int i=0; i<NMAX*Esiz; i++)
	  for (int j=i; j<NMAX*Esiz; j++)
	    SSe[0][mm][i][j] += SSe[nth][mm][i][j];

	for (int i=0; i<NMAX*Osiz; i++)
	  for (int j=i; j<NMAX*Osiz; j++)
	    SSo[0][mm][i][j] += SSo[nth][mm][i][j];
      } else {

	for (int i=0; i<NMAX*(LMAX-mm+1); i++)
	  for (int j=i; j<NMAX*(LMAX-mm+1); j++)
	    SS[0][mm][i][j] += SS[nth][mm][i][j];
      }
    }

  }

  if (VFLAG & 16) {

    for (int mm=0; mm<=MMAX; mm++) {
      bool bad = false;

      if (EvenOdd) {
	int Esiz = lE[mm].size();
	int Osiz = lO[mm].size();

	for (int i=0; i<NMAX*Esiz; i++)
	  for (int j=i; j<NMAX*Esiz; j++)
	    if (std::isnan(SCe[0][mm][i][j])) bad = true;

	for (int i=0; i<NMAX*Osiz; i++)
	  for (int j=i; j<NMAX*Osiz; j++)
	    if (std::isnan(SCo[0][mm][i][j])) bad = true;

      } else {
	for (int i=0; i<NMAX*(LMAX-mm+1); i++)
	  for (int j=i; j<NMAX*(LMAX-mm+1); j++)
	    if (std::isnan(SC[0][mm][i][j])) bad = true;
      }
      
      if (bad) {
	cerr << "Process " << myid << ": EmpCylSL.Has nan in C[" << mm << "]"
	     << endl;
      }
    }
    
    for (int mm=1; mm<=MMAX; mm++) {
      bool bad = false;

      if (EvenOdd) {
	int Esiz = lE[mm].size();
	int Osiz = lO[mm].size();

	for (int i=0; i<NMAX*Esiz; i++)
	  for (int j=i; j<NMAX*Esiz; j++)
	    if (std::isnan(SSe[0][mm][i][j])) bad = true;

	for (int i=0; i<NMAX*Osiz; i++)
	  for (int j=i; j<NMAX*Osiz; j++)
	    if (std::isnan(SSo[0][mm][i][j])) bad = true;

      } else {
	for (int i=0; i<NMAX*(LMAX-mm+1); i++)
	  for (int j=i; j<NMAX*(LMAX-mm+1); j++)
	    if (std::isnan(SS[0][mm][i][j])) bad = true;
      }
      
      if (bad) {
	cerr << "Process " << myid << ": EmpCylSL.Has nan in S[" << mm << "]"
	     << endl;
      }
    }

  }

  //
  //  Distribute covariance to all processes
  //
  for (int mm=0; mm<=MMAX; mm++) {

    if (EvenOdd) {
      int Esiz = lE[mm].size();
      int Osiz = lO[mm].size();

      icnt=0;
      for (int i=0; i<NMAX*Esiz; i++)
	for (int j=i; j<NMAX*Esiz; j++)
	  MPIin_eof[icnt++] = SCe[0][mm][i][j];
    
      MPI_Allreduce ( MPIin_eof.data(), MPIout_eof.data(), 
		      NMAX*Esiz*(NMAX*Esiz+1)/2,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      icnt=0;
      for (int i=0; i<NMAX*Esiz; i++)
	for (int j=i; j<NMAX*Esiz; j++)
	  SCe[0][mm][i][j] = MPIout_eof[icnt++];
      

      icnt=0;
      for (int i=0; i<NMAX*Osiz; i++)
	for (int j=i; j<NMAX*Osiz; j++)
	  MPIin_eof[icnt++] = SCo[0][mm][i][j];
    
      MPI_Allreduce ( MPIin_eof.data(), MPIout_eof.data(), 
		      NMAX*Osiz*(NMAX*Osiz+1)/2,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      icnt=0;
      for (int i=0; i<NMAX*Osiz; i++)
	for (int j=i; j<NMAX*Osiz; j++)
	  SCo[0][mm][i][j] = MPIout_eof[icnt++];
      
    } else {

      icnt=0;
      for (int i=0; i<NMAX*(LMAX-mm+1); i++)
	for (int j=i; j<NMAX*(LMAX-mm+1); j++)
	  MPIin_eof[icnt++] = SC[0][mm][i][j];
    
      MPI_Allreduce ( MPIin_eof.data(), MPIout_eof.data(), 
		      NMAX*(LMAX-mm+1)*(NMAX*(LMAX-mm+1)+1)/2,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      icnt=0;
      for (int i=0; i<NMAX*(LMAX-mm+1); i++)
	for (int j=i; j<NMAX*(LMAX-mm+1); j++)
	  SC[0][mm][i][j] = MPIout_eof[icnt++];
      
    }
  }
  
  for (int mm=1; mm<=MMAX; mm++) {

    if (EvenOdd) {
      int Esiz = lE[mm].size();
      int Osiz = lO[mm].size();

      icnt=0;
      for (int i=0; i<NMAX*Esiz; i++)
	for (int j=i; j<NMAX*Esiz; j++)
	  MPIin_eof[icnt++] = SSe[0][mm][i][j];
  
      MPI_Allreduce ( MPIin_eof.data(), MPIout_eof.data(), 
		      NMAX*Esiz*(NMAX*Esiz+1)/2,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      icnt=0;
      for (int i=0; i<NMAX*Esiz; i++)
	for (int j=i; j<NMAX*Esiz; j++)
	  SSe[0][mm][i][j] = MPIout_eof[icnt++];

      icnt=0;
      for (int i=0; i<NMAX*Osiz; i++)
	for (int j=i; j<NMAX*Osiz; j++)
	  MPIin_eof[icnt++] = SSo[0][mm][i][j];
  
      MPI_Allreduce ( MPIin_eof.data(), MPIout_eof.data(), 
		      NMAX*Osiz*(NMAX*Osiz+1)/2,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      icnt=0;
      for (int i=0; i<NMAX*Osiz; i++)
	for (int j=i; j<NMAX*Osiz; j++)
	  SSo[0][mm][i][j] = MPIout_eof[icnt++];

    } else {

      icnt=0;
      for (int i=0; i<NMAX*(LMAX-mm+1); i++)
	for (int j=i; j<NMAX*(LMAX-mm+1); j++)
	  MPIin_eof[icnt++] = SS[0][mm][i][j];
  
      MPI_Allreduce ( MPIin_eof.data(), MPIout_eof.data(), 
		      NMAX*(LMAX-mm+1)*(NMAX*(LMAX-mm+1)+1)/2,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      icnt=0;
      for (int i=0; i<NMAX*(LMAX-mm+1); i++)
	for (int j=i; j<NMAX*(LMAX-mm+1); j++)
	  SS[0][mm][i][j] = MPIout_eof[icnt++];
    }

  }

  //
  // DEBUG: check for nan
  //

  if (VFLAG & 16) {

    for (int n=0; n<numprocs; n++) {
      if (myid==n) {
	for (int mm=0; mm<=MMAX; mm++) {
	  bool bad = false;

	  if (EvenOdd) {
	    int Esiz = lE[mm].size();
	    int Osiz = lO[mm].size();

	    for (int i=0; i<NMAX*Esiz; i++)
	      for (int j=i; j<NMAX*Esiz; j++)
		if (std::isnan(SCe[0][mm][i][j])) bad = true;

	    for (int i=0; i<NMAX*Osiz; i++)
	      for (int j=i; j<NMAX*Osiz; j++)
		if (std::isnan(SCo[0][mm][i][j])) bad = true;

	  } else {
	    for (int i=0; i<NMAX*(LMAX-mm+1); i++)
	      for (int j=i; j<NMAX*(LMAX-mm+1); j++)
		if (std::isnan(SC[0][mm][i][j])) bad = true;
	  }
	
	  if (bad) {
	    cerr << "Process " << myid << ": EmpCylSL.Has nan in C[" << mm << "]"
		 << endl;
	  }
	}
	
	for (int mm=1; mm<=MMAX; mm++) {
	  bool bad = false;

	  if (EvenOdd) {
	    int Esiz = lE[mm].size();
	    int Osiz = lO[mm].size();

	    for (int i=0; i<NMAX*Esiz; i++)
	      for (int j=i; j<NMAX*Esiz; j++)
		if (std::isnan(SSe[0][mm][i][j])) bad = true;

	    for (int i=0; i<NMAX*Osiz; i++)
	      for (int j=i; j<NMAX*Osiz; j++)
		if (std::isnan(SSo[0][mm][i][j])) bad = true;

	  } else {
	    for (int i=0; i<NMAX*(LMAX-mm+1); i++)
	      for (int j=i; j<NMAX*(LMAX-mm+1); j++)
		if (std::isnan(SS[0][mm][i][j])) bad = true;
	  }
	  
	  if (bad) {
	    cerr << "Process " << myid << ": EmpCylSL.Has nan in S[" << mm << "]"
		 << endl;
	  }
	}
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

  }

  // END DEBUG


  if (myid==0) {

    int slave = 0;
    int request_id = 1;		// Begin with cosine case
    int M;

    M = 0;			// Initial counters
    while (M<=MMAX) {
	
      // Send request to slave
      if (slave<numprocs-1) {
	  
	slave++;
	  
	MPI_Send(&request_id, 1, MPI_INT, slave, 1, MPI_COMM_WORLD);
	MPI_Send(&M,  1, MPI_INT, slave, 2, MPI_COMM_WORLD);
	  
	// Increment counters
	request_id++;
	if (request_id>1) {
	  M++;
	  request_id = 0;
	}
	
	if (VFLAG & 16)
	  cerr << "master in make_eof: done waiting on Slave " << slave 
	       << ", next M=" << M << endl;
      }
	
				// If M>MMAX before processor queue exhausted,
				// exit loop and reap the slave data
      if (M>MMAX) break;

      if (slave == numprocs-1) {
	  
	//
	// <Wait and receive and send new request>
	//
	receive_eof(request_id, M);
	  
	// Increment counters

	request_id++;
	if (request_id>1) {
	  M++;
	  request_id = 0;
	}
      }
    }
    
    //
    // <Wait for all slaves to return and flag to continue>
    //
    if (VFLAG & 16)
      cerr << "master in make_eof: now waiting for all slaves to finish" 
	   << endl;
      
				// Dispatch resting slaves
    if (slave < numprocs-1) {
      request_id = -1;		// request_id < 0 means continue
      for (int s=numprocs-1; s>slave; s--) {
	MPI_Send(&request_id, 1, MPI_INT, s, 1, MPI_COMM_WORLD);
      }
    }

				// Get data from working slaves
    while (slave) {
      receive_eof(-1,0);
      slave--;
    }
      
  } else {

    int M, request_id;

    while (1) {
				// Wait for request . . .
      MPI_Recv(&request_id, 1, 
	       MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
      
				// Done!
      if (request_id<0) {
	if (VFLAG & 16)
	  cerr << "Slave " << setw(4) << myid 
	       << ": received DONE signal" << endl;
	break;
      }

      MPI_Recv(&M, 1, 
	       MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
	
      if (VFLAG & 16)
	cerr << "Slave " << setw(4) << myid << ": received orders type="
	     << request_id << "  M=" << M << endl;

      if (request_id) {
				// Complete symmetric part

	if (EvenOdd) {
	  int Esiz = lE[M].size();
	  int Osiz = lO[M].size();

	  for (int i=0; i<NMAX*Esiz; i++) {
	    for (int j=i; j<NMAX*Esiz; j++)
	      varE[M](i, j) = SCe[0][M][i][j];
	  }

	  for (int i=0; i<NMAX*Esiz; i++) {
	    for (int j=i; j<NMAX*Esiz; j++) {
	      varE[M](j, i) = SCe[0][M][i][j];
	    }
	  }
    
	  double maxV = 0.0;
	  for (int i=0; i<NMAX*Esiz; i++) {
	    for (int j=0; j<NMAX*Esiz; j++) {
	      tmp = fabs(varE[M](i, j));
	      if (tmp > maxV) maxV = tmp;
	    }
	  }
	  
	  varE[M] /= maxV;

	  for (int i=0; i<NMAX*Osiz; i++) {
	    for (int j=i; j<NMAX*Osiz; j++)
	      varO[M](i, j) = SCo[0][M][i][j];
	  }

	  for (int i=0; i<NMAX*Osiz; i++) {
	    for (int j=i; j<NMAX*Osiz; j++) {
	      varO[M](j, i) = SCo[0][M][i][j];
	    }
	  }
    
	  maxV = 0.0;
	  for (int i=0; i<NMAX*Osiz; i++) {
	    for (int j=i; j<NMAX*Osiz; j++) {
	      tmp = fabs(varO[M](i, j));
	      if (tmp > maxV) maxV = tmp;
	    }
	  }
	  
	  varO[M] /= maxV;

	} else {

	  for (int i=0; i<NMAX*(LMAX-M+1); i++) {
	    for (int j=i; j<NMAX*(LMAX-M+1); j++)
	      var[M](i, j) = SC[0][M][i][j];
	  }

	  for (int i=0; i<NMAX*(LMAX-M+1); i++) {
	    for (int j=i; j<NMAX*(LMAX-M+1); j++) {
	      var[M](j, i) = SC[0][M][i][j];
	    }
	  }
    
	  double maxV = 0.0;
	  for (int i=0; i<NMAX*(LMAX-M+1); i++) {
	    for (int j=i; j<NMAX*(LMAX-M+1); j++) {
	      tmp = fabs(var[M](i, j));
	      if (tmp > maxV) maxV = tmp;
	    }
	  }
	  
	  var[M] /= maxV;
	}

    
	//==========================
	// Solve eigenvalue problem 
	//==========================
    
	if (VFLAG & 16) {
	  int nancount = 0;

	  if (EvenOdd) {

	    for (int i=0; i<varE[M].rows(); i++) {
	      for (int j=0; j<varE[M].cols(); j++) {
		if (std::isnan(varE[M](i, j))) nancount++;
	      }
	    }

	    if (nancount) {
	      std::cout << "Process " << setw(4) << myid 
			<< ": in eigenvalue problem [even] with "
			<< "rank=[" << varE[M].cols() << ", " 
			<< varE[M].rows() << "]"
			<< ", found " << nancount << " NaN values" << endl;
	    }

	    nancount = 0;
	    for (int i=0; i<varO[M].rows(); i++) {
	      for (int j=0; j<varO[M].cols(); j++) {
		if (std::isnan(varO[M](i, j))) nancount++;
	      }
	    }

	    if (nancount) {
	      std::cout << "Process " << setw(4) << myid 
			<< ": in eigenvalue problem [odd] with "
			<< "rank=[" << varO[M].cols() << ", " 
			<< varO[M].rows() << "]"
			<< ", found " << nancount << " NaN values" << endl;
	    }

	  } else {

	    for (int i=0; i<var[M].rows(); i++) {
	      for (int j=0; j<var[M].cols(); j++) {
		if (std::isnan(var[M](i, j))) nancount++;
	      }
	    }
	  
	    if (nancount) {

	      std::cout << "Process " << setw(4) << myid 
			<< ": in eigenvalue problem with "
			<< "rank=[" << var[M].cols() << ", " 
			<< var[M].rows() << "]"
			<< ", found " << nancount << " NaN values" << endl;
	    }
	  }

	  timer.reset();
	  timer.start();
	}

	if (VFLAG & 32) {

	  std::ostringstream sout;
	  sout << "variance_test." << M << "." << request_id;
	  std::ofstream dout(sout.str());
	  
	  dout << std::string(60, '-') << std::endl
	       << " M=" << M           << std::endl
	       << std::string(60, '-') << std::endl;

	  if (EvenOdd) {

	    for (int i=0; i<varE[M].rows(); i++) {
	      for (int j=0; j<varE[M].cols(); j++) {
		dout << std::setw(16) << varE[M](i, j);
	      }
	      dout << std::endl;
	    }

	    dout << std::endl;
	    for (int i=0; i<varO[M].rows(); i++) {
	      for (int j=0; j<varO[M].cols(); j++) {
		dout << std::setw(16) << varO[M](i, j);
	      }
	      dout << std::endl;
	    }

	  } else {

	    for (int i=0; i<var[M].rows(); i++) {
	      for (int j=0; j<var[M].cols(); j++) {
		dout << std::setw(16) << var[M](i, j);
	      }
	      dout << std::endl;
	    }
	  }
	}

	if (EvenOdd) {

	  Eigen::EigenSolver<Eigen::MatrixXd> esE(varE[M]);
	  Eigen::EigenSolver<Eigen::MatrixXd> esO(varO[M]);

	  evE = esE.eigenvalues().real();
	  evO = esO.eigenvalues().real();

	  efE = esE.eigenvectors().real();
	  efO = esO.eigenvectors().real();

	  if (VFLAG & 32) {

	    std::ostringstream sout;
	    sout << "ev_test." << M << "." << request_id;
	    std::ofstream dout(sout.str());
	  
	    for (int i=0; i<efE.cols(); i++) {
	      for (int j=0; j<efE.cols(); j++) {
		double sum = efE.col(i).adjoint() * efE.col(j);
		dout << std::setw(4) << i << std::setw(4) << j
		     << std::setw(18) << sum << std::endl;
	      }
	    }

	    dout << std::endl;
	  
	    for (int i=0; i<efO.cols(); i++) {
	      for (int j=0; j<efO.cols(); j++) {
		double sum = efO.col(i).adjoint() * efO.col(j);
		dout << std::setw(4) << i << std::setw(4) << j
		     << std::setw(18) << sum << std::endl;
	      }
	    }
	  }

	} else {

	  Eigen::EigenSolver<Eigen::MatrixXd> es(var[M]);

	  ev = es.eigenvalues().real();
	  ef = es.eigenvectors().real();

	  if (VFLAG & 32) {

	    std::ostringstream sout;
	    sout << "ev_test." << M << "." << request_id;
	    std::ofstream dout(sout.str());
	  
	    for (int i=0; i<ef.cols(); i++) {
	      for (int j=0; j<ef.cols(); j++) {
		double sum = ef.col(i).adjoint() * ef.col(j);
		dout << std::setw(4) << i << std::setw(4) << j
		     << std::setw(18) << sum << std::endl;
	      }
	    }
	  }
	}

	if (VFLAG & 16) {
	  cout << "Process " << setw(4) << myid 
	       << ": completed eigenproblem in " 
	       << timer.stop() << " seconds"
	       << endl;
	}

      } else {
				// Complete symmetric part
    
	if (EvenOdd) {
	  int Esiz = lE[M].size();
	  int Osiz = lO[M].size();

	  for (int i=0; i<NMAX*Esiz; i++) {
	    for (int j=i; j<NMAX*Esiz; j++)
	      varE[M](i, j) = SSe[0][M][i][j];
	  }

	  for (int i=0; i<NMAX*Esiz; i++) {
	    for (int j=i; j<NMAX*Esiz; j++) {
	      varE[M](j, i) = SSe[0][M][i][j];
	    }
	  }
    
	  double maxV = 0.0;
	  for (int i=0; i<NMAX*Esiz; i++) {
	    for (int j=i; j<NMAX*Esiz; j++) {
	      tmp = fabs(varE[M](i, j));
	      if (tmp > maxV) maxV = tmp;
	    }
	  }
	  
	  if (maxV>1.0e-5)
	    varE[M] /= maxV;
    
	  for (int i=0; i<NMAX*Osiz; i++) {
	    for (int j=i; j<NMAX*Osiz; j++)
	      varO[M](i, j) = SSo[0][M][i][j];
	  }

	  for (int i=0; i<NMAX*Osiz; i++) {
	    for (int j=i; j<NMAX*Osiz; j++) {
	      varO[M](j, i) = SSo[0][M][i][j];
	    }
	  }
    
	  maxV = 0.0;
	  for (int i=0; i<NMAX*Osiz; i++) {
	    for (int j=i; j<NMAX*Osiz; j++) {
	      tmp = fabs(varO[M](i, j));
	      if (tmp > maxV) maxV = tmp;
	    }
	  }
	  
	  if (maxV>1.0e-5)
	    varO[M] /= maxV;
    
	} else {

	  for (int i=0; i<NMAX*(LMAX-M+1); i++) {
	    for (int j=i; j<NMAX*(LMAX-M+1); j++)
	      var[M](i, j) = SS[0][M][i][j];
	  }
	  
	  for (int i=0; i<NMAX*(LMAX-M+1); i++) {
	    for (int j=i; j<NMAX*(LMAX-M+1); j++) {
	      var[M](j, i) = SS[0][M][i][j];
	    }
	  }
	  
	  double maxV = 0.0;
	  for (int i=0; i<NMAX*(LMAX-M+1); i++) {
	    for (int j=i; j<NMAX*(LMAX-M+1); j++) {
	      tmp = fabs(var[M](i, j));
	      if (tmp > maxV) maxV = tmp;
	    }
	  }
	  
	  if (maxV>1.0e-5)
	    var[M] /= maxV;
	}

	//==========================
	// Solve eigenvalue problem
	//==========================
    
	if (VFLAG & 16) {

	  int nancount = 0;

	  if (EvenOdd) {

	    for (int i=0; i<varE[M].rows(); i++) {
	      for (int j=0; j<varE[M].cols(); j++) {
		if (std::isnan(varE[M](i, j))) nancount++;
	      }
	    }
	    
	    cout << "Process " << setw(4) << myid 
		 << ": in eigenvalue problem [even] with "
		 << "rank=[" << varE[M].cols() << ", " 
		 << varE[M].rows() << "]"
		 << ", found " << nancount << " NaN values" << endl;

	    nancount = 0;
	    for (int i=0; i<varO[M].cols(); i++) {
	      for (int j=0; j<varO[M].cols(); j++) {
		if (std::isnan(varO[M](i, j))) nancount++;
	      }
	    }
	    
	    cout << "Process " << setw(4) << myid 
		 << ": in eigenvalue problem [odd] with "
		 << "rank=[" << varO[M].cols() << ", " 
		 << varO[M].rows() << "]"
		 << ", found " << nancount << " NaN values" << endl;

	  } else {

	    for (int i=0; i<var[M].rows(); i++) {
	      for (int j=0; j<var[M].cols(); j++) {
		if (std::isnan(var[M](i, j))) nancount++;
	      }
	    }
	    
	    cout << "Process " << setw(4) << myid 
		 << ": in eigenvalue problem with "
		 << "rank=[" << var[M].cols() << ", " 
		 << var[M].rows() << "]"
		 << ", found " << nancount << " NaN values" << endl;
	  }

	  timer.reset();
	  timer.start();
	}

	if (EvenOdd) {

	  Eigen::EigenSolver<Eigen::MatrixXd> esE(varE[M]);
	  Eigen::EigenSolver<Eigen::MatrixXd> esO(varO[M]);

	  evE = esE.eigenvalues().real();
	  evO = esO.eigenvalues().real();

	  efE = esE.eigenvectors().real();
	  efO = esO.eigenvectors().real();

	  if (VFLAG & 32) {

	    std::ostringstream sout;
	    sout << "ev_test." << M << "." << request_id;
	    std::ofstream dout(sout.str());
	  
	    for (int i=0; i<efE.cols(); i++) {
	      for (int j=0; j<efE.cols(); j++) {
		double sum = efE.col(i).adjoint() * efE.col(j);
		dout << std::setw(4) << i << std::setw(4) << j
		     << std::setw(18) << sum << std::endl;
	      }
	    }

	    dout << std::endl;
	  
	    for (int i=0; i<efO.cols(); i++) {
	      for (int j=0; j<efO.cols(); j++) {
		double sum = efO.col(i).adjoint() * efO.col(j);
		dout << std::setw(4) << i << std::setw(4) << j
		     << std::setw(18) << sum << std::endl;
	      }
	    }
	  }

	} else {

	  Eigen::EigenSolver<Eigen::MatrixXd> es(var[M]);

	  ev = es.eigenvalues().real();
	  ef = es.eigenvectors().real();

	  if (VFLAG & 32) {

	    std::ostringstream sout;
	    sout << "ev_test." << M << "." << request_id;
	    std::ofstream dout(sout.str());
	  
	    for (int i=0; i<ef.cols(); i++) {
	      for (int j=0; j<ef.cols(); j++) {
		double sum = ef.col(i).adjoint() * ef.col(j);
		dout << std::setw(4) << i << std::setw(4) << j
		     << std::setw(18) << sum << std::endl;
	      }
	    }
	  }
	}

	if (VFLAG & 16) {
	  cout << "Process " << setw(4) << myid 
	       << ": completed eigenproblem in " 
	       << timer.stop() << " seconds"
	       << endl;
	}

      }

      if (VFLAG & 2)
	cerr << "Slave " << setw(4) << myid 
	     << ": with request_id=" << request_id
	     << ", M=" << M << " calling compute_eof_grid" << endl;

      if (VFLAG & 16) {
	timer.reset();
	timer.start();
      }

      if (EvenOdd)
	compute_even_odd(request_id, M);
      else
	compute_eof_grid(request_id, M);

      if (VFLAG & 16) {
	cout << "Process " << setw(4) << myid << ": completed EOF grid for id="
	     << request_id << " and M=" << M << " in " 
	     << timer.stop() << " seconds"
	     << endl;
      }
      else if (VFLAG & 2)
	cerr << "Slave " << setw(4) << myid 
	     << ": with request_id=" << request_id
	     << ", M=" << M << " COMPLETED compute_eof_grid" << endl;
    }

  }
				// Send grid to all processes
  if (VFLAG & 2) {
    MPI_Barrier(MPI_COMM_WORLD);
    cerr << "Process " << setw(4) << myid 
	 << ": about to enter send_eof_grid" << endl;
  }

  if (VFLAG & 16) {
    timer.reset();
    timer.start();
  }

  send_eof_grid();

  if (VFLAG & 16) {
    cout << "Process " << setw(4) << myid << ": grid reduced in " 
	 << timer.stop()  << " seconds"
	 << endl;
  } 
  else if (VFLAG & 2)
    cerr << "Process " << setw(4) << myid << ": grid reduce completed" << endl;

				// Cache table for restarts
				// (it would be nice to multithread or fork
				//  this call . . . )
  if (myid==0) cache_grid(1);
  
  eof_made      = true;
  coefs_made    = std::vector<short>(multistep+1, false);

  if (VFLAG & 2) {
    MPI_Barrier(MPI_COMM_WORLD);
    cerr << "Process " << setw(4) << myid 
	 << ": EOF computation completed" << endl;
  }
}


void EmpCylSL::accumulate_eof(std::vector<Particle>& part, bool verbose)
{

  double r, phi, z, mass;

  int ncnt=0;
  if (myid==0 && verbose) cout << endl;

  setup_eof();

  for (auto p=part.begin(); p!=part.end(); p++) {

    mass = p->mass;
    r = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
    phi = atan2(p->pos[1], p->pos[0]);
    z = p->pos[2];
    
    accumulate_eof(r, z, phi, mass, 0, p->level);
    if (myid==0 && verbose) {
      if ( (ncnt % 100) == 0) cout << "\r>> " << ncnt << " <<" << flush;
      ncnt++;
    }
  }

}
  

void EmpCylSL::accumulate_eof_thread(std::vector<Particle>& part, bool verbose)
{
  setup_eof();

  std::vector<std::thread> t(nthrds);
 
  // Launch the threads
  for (int id=0; id<nthrds; ++id) {
    t[id] = std::thread(&EmpCylSL::accumulate_eof_thread_call, this, id, &part, verbose);
  }
  // Join the threads
  for (int id=0; id<nthrds; ++id) {
    t[id].join();
  }
}


void EmpCylSL::accumulate_eof_thread_call(int id, std::vector<Particle>* p, bool verbose)
{
  int nbodies = p->size();
    
  if (nbodies == 0) return;

  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double r, phi, z, mass;

  int ncnt=0;
  if (myid==0 && id==0 && verbose) cout << endl;
  
  for (int n=nbeg; n<nend; n++) {
				// Phase space coords
    mass = (*p)[n].mass;
    r    = sqrt((*p)[n].pos[0]*(*p)[n].pos[0] + (*p)[n].pos[1]*(*p)[n].pos[1]);
    phi  = atan2((*p)[n].pos[1], (*p)[n].pos[0]);
    z    = (*p)[n].pos[2];
				// Call accumulation for this particle
    accumulate_eof(r, z, phi, mass, id, (*p)[n].level);

    if (myid==0 && id==0 && verbose) {
      if ( (ncnt % 100) == 0) cout << "\r>> " << ncnt << " <<" << flush;
      ncnt++;
    }
  }

}
  

void EmpCylSL::accumulate(std::vector<Particle>& part, int mlevel,
			  bool verbose, bool compute)
{
   double r, phi, z, mass;

  int ncnt=0;
  if (myid==0 && verbose) cout << endl;

  setup_accumulation();

  for (auto p=part.begin(); p!=part.end(); p++) {

    double mass = p->mass;
    double r    = sqrt(p->pos[0]*p->pos[0] + p->pos[1]*p->pos[1]);
    double phi  = atan2(p->pos[1], p->pos[0]);
    double z    = p->pos[2];
    
    accumulate(r, z, phi, mass, p->indx, 0, mlevel, compute);

    if (myid==0 && verbose) {
      if ( (ncnt % 100) == 0) cout << "\r>> " << ncnt << " <<" << flush;
      ncnt++;
    }
  }

}
  

void EmpCylSL::accumulate_thread(std::vector<Particle>& part, int mlevel, bool verbose)
{
  setup_accumulation();

  std::vector<std::thread> t(nthrds);
 
  // Launch the threads
  //
  for (int id=0; id<nthrds; ++id) {
    t[id] = std::thread(&EmpCylSL::accumulate_thread_call, this, id, &part, mlevel, verbose);
  }

  // Join the threads
  //
  for (int id=0; id<nthrds; ++id) {
    t[id].join();
  }
}


void EmpCylSL::accumulate_thread_call(int id, std::vector<Particle>* p, int mlevel, bool verbose)
{
  int nbodies = p->size();
    
  if (nbodies == 0) return;

  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  int ncnt=0;
  if (myid==0 && id==0 && verbose) cout << endl;

  for (int n=nbeg; n<nend; n++) {
    
    double mass = (*p)[n].mass;
    double r    = sqrt((*p)[n].pos[0]*(*p)[n].pos[0] + (*p)[n].pos[1]*(*p)[n].pos[1]);
    double phi  = atan2((*p)[n].pos[1], (*p)[n].pos[0]);
    double z    = (*p)[n].pos[2];
    
    accumulate(r, z, phi, mass, (*p)[n].indx, id, mlevel);

    if (myid==0 && id==0 && verbose) {
      if ( (ncnt % 100) == 0) cout << "\r>> " << ncnt << " <<" << flush;
      ncnt++;
    }
  }

}
  

void EmpCylSL::getPotParticle(double x, double y, double z,
			      EmpCylSL::ContribArray &retC,
			      EmpCylSL::ContribArray &retS)
{
  constexpr double norm = -4.0*M_PI;

  if (retC.size()==0) {
    retC.resize(MMAX+1);
    retS.resize(MMAX+1);
    for (auto & v : retC) v.resize(rank3);
    for (auto & v : retS) v.resize(rank3);
  }

  double R = sqrt(x*x + y*y);
  double phi = atan2(y, x);
  
  get_pot(vc[0], vs[0], R, z);
  for (int mm=0; mm<=MMAX; mm++) {
    for (int nn=0; nn<rank3; nn++) {
      retC[mm][nn] = vc[0](mm, nn)*cos(phi*mm);
      if (mm>0) retS[mm][nn] = vs[0](mm, nn)*sin(phi*mm);
    }
  }
}


void EmpCylSL::accumulate(double r, double z, double phi, double mass, 
			  unsigned long seq, int id, int mlevel, bool compute)
{

  if (coefs_made[mlevel]) {
    ostringstream ostr;
    ostr << "EmpCylSL::accumulate: Process " << myid << ", Thread " << id 
	 << ": calling setup_accumulation from accumulate, aborting" << endl;
    throw GenericError(ostr.str(), __FILE__, __LINE__);
  }

  double rr = sqrt(r*r+z*z);
  if (rr/ASCALE>Rtable) return;

  howmany1[mlevel][id]++;

  double msin, mcos;
  int mm;
  
  double norm = -4.0*M_PI;
  
  unsigned whch;
  if (compute and PCAVAR) {
    whch = seq % sampT;
    pthread_mutex_lock(&used_lock);
    numbT1[id][whch] += 1;
    massT1[id][whch] += mass;
    pthread_mutex_unlock(&used_lock);
  }

  get_pot(vc[id], vs[id], r, z);

  for (mm=0; mm<=MMAX; mm++) {

    mcos = cos(phi*mm);
    msin = sin(phi*mm);

    for (int nn=0; nn<rank3; nn++) {
      double hold = norm * mass * mcos * vc[id](mm, nn);

      cosN(mlevel)[id][mm][nn] += hold;

      if (compute and PCAVAR) {
	double hc1 = vc[id](mm, nn)*mcos, hs1 = 0.0;
	if (mm) hs1 = vs[id](mm, nn)*msin;
	double modu1 = sqrt(hc1*hc1 + hs1*hs1);

	covV(id, whch, mm)[nn] += mass * modu1;

	for (int oo=0; oo<rank3; oo++) {
	  double hc2 = vc[id](mm, oo)*mcos, hs2 = 0.0;
	  if (mm) hs2 = vs[id](mm, oo)*msin;
	  double modu2 = sqrt(hc2*hc2 + hs2*hs2);

	  covM(id, whch, mm)(nn, oo) += mass * modu1 * modu2;
	}
      }

      if (mm>0) {
	hold = norm * mass * msin * vs[id](mm, nn);
	sinN(mlevel)[id][mm][nn] += hold;
      }

      if (compute and PCAEOF) {
	double hold1 = vc[id](mm, nn), hold2 = 0.0;
	if (mm>0) hold2 = vs[id](mm, nn);
	double modu1 = sqrt(hold1*hold1 + hold2*hold2)*norm;
	for (int oo=0; oo<rank3; oo++) {
	  hold1 = vc[id](mm, oo), hold2 = 0.0;
	  if (mm>0) hold2 = vs[id](mm, oo);
	  double modu2 = sqrt(hold1*hold1 + hold2*hold2)*norm;
	  tvar[id][mm](nn, oo) += modu1 * modu2 * mass;
	}
      }
    }

    cylmass1[id] += mass;
  }

}


void EmpCylSL::make_coefficients(unsigned M0, bool compute)
{
  if (MPIin.size()==0) {
				// Vector reduction
    MPIin  .resize(rank3*(MMAX+1));
    MPIout .resize(rank3*(MMAX+1));
				// Matrix reduction
    MPIin2 .resize(rank3*rank3*(MMAX+1));
    MPIout2.resize(rank3*rank3*(MMAX+1));
  }
  

  for (unsigned M=M0; M<=multistep; M++) {
    
    if (coefs_made[M]) continue;

				// Sum up over threads
				//
    for (int nth=1; nth<nthrds; nth++) {

      howmany1[M][0] += howmany1[M][nth];

      for (int mm=0; mm<=MMAX; mm++)
	for (int nn=0; nn<rank3; nn++) {
	  cosN(M)[0][mm][nn] += cosN(M)[nth][mm][nn];
	}
      
      for (int mm=1; mm<=MMAX; mm++)
	for (int nn=0; nn<rank3; nn++) {
	  sinN(M)[0][mm][nn] += sinN(M)[nth][mm][nn];
	}
    }
				// Begin distribution loop
				//
    for (int mm=0; mm<=MMAX; mm++)
      for (int nn=0; nn<rank3; nn++)
	MPIin[mm*rank3 + nn] = cosN(M)[0][mm][nn];
    
    MPI_Allreduce ( MPIin.data(), MPIout.data(), rank3*(MMAX+1),
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int mm=0; mm<=MMAX; mm++)
      for (int nn=0; nn<rank3; nn++)
	if (multistep)
	  cosN(M)[0][mm][nn] = MPIout[mm*rank3 + nn];
	else
	  accum_cos[mm][nn] = MPIout[mm*rank3 + nn];
    

    for (int mm=1; mm<=MMAX; mm++)
      for (int nn=0; nn<rank3; nn++)
	MPIin[mm*rank3 + nn] = sinN(M)[0][mm][nn];
    
    MPI_Allreduce ( MPIin.data(), MPIout.data(), rank3*(MMAX+1),
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

    for (int mm=1; mm<=MMAX; mm++)
      for (int nn=0; nn<rank3; nn++)
	if (multistep)
	  sinN(M)[0][mm][nn] = MPIout[mm*rank3 + nn];
	else
	  accum_sin[mm][nn] = MPIout[mm*rank3 + nn];
    
    coefs_made[M] = true;
  }
  

  if (compute) {
				// Sum up over threads
				//
    for (int nth=1; nth<nthrds; nth++) {

      if (PCAEOF) {
	for (int mm=0; mm<=MMAX; mm++)
	  for (int nn=0; nn<rank3; nn++)
	    for (int oo=0; oo<rank3; oo++)
	      tvar[0][mm](nn, oo) += tvar[nth][mm](nn, oo);
      }
	
      for (unsigned T=0; T<sampT; T++) {
	numbT1[0][T] += numbT1[nth][T];
	massT1[0][T] += massT1[nth][T];

	for (int mm=0; mm<=MMAX; mm++) {
	  for (int nn=0; nn<rank3; nn++) {
	    covV(0, T, mm)[nn] += covV(nth, T, mm)[nn];
	    for (int oo=0; oo<rank3; oo++) {
	      covM(0, T, mm)(nn, oo) += covM(nth, T, mm)(nn, oo);
	    }
	  }
	}

      } // T loop
      
    } // Thread loop


    // Mass used to compute variance in each partition
    //
    if (PCAVAR) {
      MPI_Allreduce ( &numbT1[0][0], &numbT[0], sampT,
		      MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce ( &massT1[0][0], &massT[0], sampT,
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    // Test variance
    //
    if (PCAEOF) {

      std::vector<double> MPIinT(rank3*rank3*(1 + MMAX));
      std::vector<double> MPIotT(rank3*rank3*(1 + MMAX));
    
      for (int mm=0; mm<=MMAX; mm++)
	for (int nn=0; nn<rank3; nn++)
	  for (int oo=0; oo<rank3; oo++)
	    MPIinT[mm*rank3*rank3 + nn*rank3 + oo] = tvar[0][mm](nn, oo);
  
      MPI_Allreduce ( MPIinT.data(), MPIotT.data(), rank3*rank3*(MMAX+1),
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      
      for (int mm=0; mm<=MMAX; mm++)
	for (int nn=0; nn<rank3; nn++)
	  for (int oo=0; oo<rank3; oo++)
	    tvar[0][mm](nn, oo) = MPIotT[mm*rank3*rank3 + nn*rank3 + oo];
    }
      
    // Begin distribution loop for variance jackknife
    //
    for (unsigned T=0; T<sampT; T++) {
      
      for (int mm=0; mm<=MMAX; mm++) {
	for (int nn=0; nn<rank3; nn++) {
	  MPIin[mm*rank3 + nn] = covV(0, T, mm)[nn];
	  for (int oo=0; oo<rank3; oo++) {
	    MPIin2[mm*rank3*rank3 + nn*rank3 + oo] = covM(0, T, mm)(nn, oo);
	  }
	}
      }

      MPI_Allreduce ( MPIin.data(), MPIout.data(), rank3*(MMAX+1),
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce ( MPIin2.data(), MPIout2.data(), rank3*rank3*(MMAX+1),
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for (int mm=0; mm<=MMAX; mm++) {
	for (int nn=0; nn<rank3; nn++) {
	  covV(0, T, mm)[nn] = MPIout[mm*rank3 + nn];
	  for (int oo=0; oo<rank3; oo++) {
	    covM(0, T, mm)(nn, oo) = MPIout2[mm*rank3*rank3 + nn*rank3 + oo];
	  }
	}
      }
    }
    // END: T loop
  }
}

void EmpCylSL::reset_mass(void)
{ 
  cylmass=0.0; 
  cylmass_made=false; 
  for (int n=0; n<nthrds; n++) cylmass1[n] = 0.0;
}

void EmpCylSL::make_coefficients(bool compute)
{
  if (!cylmass_made) {
    MPI_Allreduce(&cylmass1[0], &cylmass, 1, MPI_DOUBLE, MPI_SUM, 
		  MPI_COMM_WORLD);
    cylmass_made = true;
  }
  
  if (coefs_made_all()) return;

  if (MPIin.size()==0) {	// Vector reduction
    MPIin  .resize(rank3*(MMAX+1));
    MPIout .resize(rank3*(MMAX+1));
				// Matrix reduction
    MPIin2 .resize(rank3*rank3*(MMAX+1));
    MPIout2.resize(rank3*rank3*(MMAX+1));
  }
				// Sum up over threads
				// 
  for (unsigned M=0; M<=multistep; M++) {

    if (coefs_made[M]) continue;
    
    for (int nth=1; nth<nthrds; nth++) {

      howmany1[M][0] += howmany1[M][nth];

      if (compute and PCAVAR) {
	for (unsigned T=0; T<sampT; T++) {
	  numbT1[0][T] += numbT1[nth][T];
	  massT1[0][T] += massT1[nth][T];
	}
      }
      
      for (int mm=0; mm<=MMAX; mm++) {
	for (int nn=0; nn<rank3; nn++) {
	  cosN(M)[0][mm][nn] += cosN(M)[nth][mm][nn];
	  if (compute and PCAVAR) {
	    for (unsigned T=0; T<sampT; T++) {
	      covV(0, T, mm)[nn] += covV(nth, T, mm)[nn];
	      for (int oo=0; oo<rank3; oo++) {
		covM(0, T, mm)(nn, oo) += covM(nth, T, mm)(nn, oo);
	      }
	    }
	  }
	  // END: PCAVAR
	}
      }

      for (int mm=1; mm<=MMAX; mm++) {
	for (int nn=0; nn<rank3; nn++) {
	  sinN(M)[0][mm][nn] += sinN(M)[nth][mm][nn];
	}
      }

    }
  }

				// Begin distribution loop

  howmany. resize(multistep+1);	// Resize and zero, if empty
  howmany1.resize(multistep+1);
  for (unsigned M=0; M<=multistep; M++) {
    for (int nth=0; nth<nthrds; nth++) howmany1[M].resize(nthrds, 0);
  }

  for (unsigned M=0; M<=multistep; M++) {

    if (coefs_made[M]) continue;

				// "howmany" is only used for debugging
    MPI_Allreduce ( &howmany1[M][0], &howmany[M], 1, MPI_UNSIGNED,
		    MPI_SUM, MPI_COMM_WORLD);

    for (int mm=0; mm<=MMAX; mm++)
      for (int nn=0; nn<rank3; nn++)
	MPIin[mm*rank3 + nn] = cosN(M)[0][mm][nn];
  
    MPI_Allreduce ( MPIin.data(), MPIout.data(), rank3*(MMAX+1),
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (int mm=0; mm<=MMAX; mm++)
      for (int nn=0; nn<rank3; nn++)
	if (multistep)
	  cosN(M)[0][mm][nn] = MPIout[mm*rank3 + nn];
	else
	  accum_cos[mm][nn] = MPIout[mm*rank3 + nn];
  }
  

  if (compute and PCAVAR) {

    for (unsigned T=0; T<sampT; T++) {
      for (int mm=0; mm<=MMAX; mm++) {
	for (int nn=0; nn<rank3; nn++) {
	  MPIin[mm*rank3 + nn] = covV(0, T, mm)[nn];
	  for (int oo=0; oo<rank3; oo++) {
	    MPIin2[mm*rank3*rank3 + nn*rank3 + oo] = covM(0, T, mm)(nn, oo);
	  }
	}
      }

      MPI_Allreduce ( MPIin.data(), MPIout.data(), rank3*(MMAX+1),
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      MPI_Allreduce ( MPIin2.data(), MPIout2.data(), rank3*rank3*(MMAX+1),
		      MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for (int mm=0; mm<=MMAX; mm++) {
	for (int nn=0; nn<rank3; nn++) {
	  covV(0, T, mm)[nn] = MPIout[mm*rank3 + nn];
	  for (int oo=0; oo<rank3; oo++) {
	    covM(0, T, mm)(nn, oo) = MPIout2[mm*rank3*rank3 + nn*rank3 + oo];
	  }
	}
      }

    } // T loop

  } // SELECT


  for (unsigned M=0; M<=multistep; M++) {
    
    if (coefs_made[M]) continue;

    for (int mm=1; mm<=MMAX; mm++)
      for (int nn=0; nn<rank3; nn++)
	MPIin[mm*rank3 + nn] = sinN(M)[0][mm][nn];
  
    MPI_Allreduce ( MPIin.data(), MPIout.data(), rank3*(MMAX+1),
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  

    for (int mm=1; mm<=MMAX; mm++)
      for (int nn=0; nn<rank3; nn++)
	if (multistep)
	  sinN(M)[0][mm][nn] = MPIout[mm*rank3 + nn];
	else
	  accum_sin[mm][nn] = MPIout[mm*rank3 + nn];
  }
  
  if (compute and PCAVAR) {
				// Mass used to compute variance in
				// each partition

    MPI_Allreduce ( &numbT1[0][0], &numbT[0], sampT,
		    MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce ( &massT1[0][0], &massT[0], sampT,
		    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
  } // END: 'compute' stanza
  
  coefs_made = vector<short>(multistep+1, true);

}


void EmpCylSL::pca_hall(bool compute, bool subsamp)
{
  if (VFLAG & 4)
    cerr << "Process " << setw(4) << myid << ": made it to pca_hall" << endl;
  
  if (compute) {

    // For PCA jack knife
    //
    if (pb)
      pb->reset();
    else
      pb = PCAbasisPtr(new PCAbasis(MMAX, rank3));
    
    static unsigned ocount = 0, eofcount = 0;

#if defined HAVE_VTK
    VtkPCAptr vtkpca;

    if (PCAVAR and PCAVTK and myid==0) {

      if (ocount==0) {	       // Look for restart position; this is
	while (1) {	       // time consuming but is only done once.
	  std::ostringstream fileN;
	  fileN << hallfile << "_"
		<< std::setfill('0') << std::setw(5) << ocount << ".vtr";
	  std::cout << "File: " << fileN.str() << std::endl;
	  std::ifstream infile(fileN.str());
	  if (not infile.good()) break;
	  ocount++;
	}
	if (ocount)
	  std::cout << "Restart in EmpCylSL::pca_hall: "
		    << "vtk output will begin at "
		    << ocount << std::endl;
      }
      
      if (ocount % VTKFRQ==0) vtkpca = VtkPCAptr(new VtkPCA(rank3));
    }
#endif

    if (PCAVAR) {
      for (auto v : numbT) pb->Tnumb += v;
      for (auto v : massT) pb->Tmass += v;
    }
    
    if (VFLAG & 4)
      cerr << "Process " << setw(4) << myid << ": mass " 
	   << pb->Tmass << " of particles" << endl;
    
    // No data?
    //
    if (PCAVAR) {
      if (pb->Tmass<=0.0) return;	
    }
    
    // Setup for diagnostic output
    //
    std::ofstream hout, mout;
    const int fwid = 14;
    if (myid==0 && hallfile.length()>0) {
      std::ostringstream ofile, mfile;
      ofile << hallfile << ".pcalog";
      hout.open(ofile.str(), ios::out | ios::app);
      if (hout.good()) {
	hout << "#" << endl << std::right
	     << "#" << endl
	     << setw( 4) << "m" << setw(4) << "n"
	     << setw(fwid) << "coef";
	if (PCAVAR) {
	  hout << setw(fwid) << "|coef|^2"
	       << setw(fwid) << "var(coef)"
	       << setw(fwid) << "cum var"
	       << setw(fwid) << "S/N (1p)"
	       << setw(fwid) << "b_Hall"
	       << setw(fwid) << "s_Hall"
	       << setw(fwid) << "|coef|_0^2"
	       << setw(fwid) << "var_0"
	       << setw(fwid) << "S/N (Np)";
	}
	if (PCAEOF) hout << setw(fwid) << "EOF";
	hout << std::endl << std::endl;
      } else {
	cerr << "Could not open <" << ofile.str() << "> for appending output" 
	     << endl;
      }
    
      mfile << hallfile << ".pcamat";
      mout.open(mfile.str(), ios::out | ios::app);
      if (mout.good()) {
	mout << "#" << endl << std::right
	     << "#" << endl << setprecision(4);
      } else {
	cerr << "Could not open <" << mfile.str() << "> for appending output" 
	     << endl;
      }
    }
    
    Eigen::VectorXd eofvec;

    if (PCAEOF) eofvec.resize(rank3);

    // Loop through each harmonic subspace [EVEN cosines]
    //
    for (int mm=0; mm<=MMAX; mm++) {
      
      if (PCAVAR) {

	// Data partitions for variance
	//
	for (unsigned T=0; T<sampT; T++) {
	  
	  if (massT[T] <= 0.0) continue; // Skip empty partition
	
	  for (int nn=0; nn<rank3; nn++) {

	    // Compute mean coef from subsample
	    //
	    (*pb)[mm]->meanJK[nn] += covV(0, T, mm)[nn] / massT[T] / sampT;
	    
	    if (subsamp) {

	      for (int oo=0; oo<rank3; oo++) {

		// Compute mean squared coefs from subsample
		//
		(*pb)[mm]->covrJK(nn, oo) +=
		  covV(0, T, mm)[nn] / massT[T] *
		  covV(0, T, mm)[oo] / massT[T] / sampT;
	      }
	      
	    } else {

	      for (int oo=0; oo<rank3; oo++) {

		// Compute mean squared coefs from subsample
		//
		(*pb)[mm]->covrJK(nn, oo) += covM(0, T, mm)(nn, oo) / massT[T] / sampT;
	      }
	    }
	  }
	}
	
	// Compute subsample variance
	//
	for (int nn=0; nn<rank3; nn++) { 
	  for (int oo=0; oo<rank3; oo++) {
	    (*pb)[mm]->covrJK(nn, oo) -=
	      (*pb)[mm]->meanJK[nn] * (*pb)[mm]->meanJK[oo];
	  }
	}
	
	Eigen::EigenSolver<Eigen::MatrixXd> es((*pb)[mm]->covrJK);

	(*pb)[mm]->evalJK = es.eigenvalues().real();
	(*pb)[mm]->evecJK = es.eigenvectors().real();
      }
    
      // Transformation output
      //
      if (mout.good()) {
	mout << "#" << std::endl
	     << "# m=" << mm << std::endl
	     << "#" << std::endl;
	if (PCAVAR) {
	  mout << "# Eigenvalues"  << std::endl;
	  double enorm = 0.0, ecum = 0.0;
	  for (int nn=0; nn<rank3; nn++)
	    enorm += (*pb)[mm]->evalJK[nn];
	  for (int nn=0; nn<rank3; nn++) {
	    ecum += (*pb)[mm]->evalJK[nn];
	    mout << std::setw( 4) << nn
		 << std::setw(12) << (*pb)[mm]->evalJK[nn]
		 << std::setw(12) << ecum / enorm
		 << std::endl;
	  }
	  mout << "# Eigenvectors" << std::endl;
	  for (int nn=0; nn<rank3; nn++) {
	    for (int oo=0; oo<rank3; oo++) {
	      mout << std::setw(12) << (*pb)[mm]->evecJK.col(nn)[oo];
	    }
	    mout << std::endl;
	  }
	  
	  mout << "# Norms" << std::endl;
	  for (int nn=0; nn<rank3; nn++) {
	    for (int oo=0; oo<rank3; oo++) {
	      mout << std::setw(12) <<
		(*pb)[mm]->evecJK.row(nn).adjoint() *
		(*pb)[mm]->evecJK.row(oo);

	    }
	    mout << std::endl;
	  }

	  mout << "# Covariance matrix" << std::endl;
	  for (int nn=0; nn<rank3; nn++) {
	    for (int oo=0; oo<rank3; oo++)
	      mout << std::setw(12) << (*pb)[mm]->covrJK(nn, oo);
	    mout << std::endl;
	  }
	}

	if (PCAEOF) {
	  Eigen::VectorXd evalVar(rank3);
	  Eigen::MatrixXd evecVar(rank3, rank3);
	  Eigen::EigenSolver<Eigen::MatrixXd> es(tvar[0][mm]);

	  evalVar = es.eigenvalues().real();
	  evecVar = es.eigenvectors().real();

	  mout << "# EOF eigenvalues" << std::endl;
	  double total = 0.0;
	  for (int nn=0; nn<rank3; nn++) {
	    total += evalVar[nn];
	    mout << std::setw(12) << evalVar[nn];
	  }
	  mout << std::endl;
	  
	  mout << "# EOF accumulation" << std::endl;
	  double cum = 0.0;
	  for (int nn=0; nn<rank3; nn++) {
	    cum += evalVar[nn];
	    mout << std::setw(12) << cum/total;
	  }
	  mout << std::endl;
	  
	  mout << "# EOF eigenvectors" << std::endl;
	  for (int nn=0; nn<rank3; nn++) {
	    for (int oo=0; oo<rank3; oo++)
	      mout << std::setw(12) << evecVar.col(nn)[oo];
	    mout << std::endl;
	  }
	  
	  Eigen::VectorXd initVar(rank3);
	  for (int nn=0; nn<rank3; nn++) {
	    initVar[nn] = accum_cos[mm][nn] * accum_cos[mm][nn];
	    if (mm) initVar[nn] += accum_sin[mm][nn] * accum_sin[mm][nn];
	    initVar[nn] = sqrt(initVar[nn]);
	  }

	  eofvec = evecVar.transpose() * initVar;

	  // VTK basis
	  //
	  for (int nn=0; nn<rank3; nn++) {
	    dump_images_basis_eof(runtag, 0.1, 0.01, 100, 40, mm, nn, eofcount,
				  evecVar.col(nn));
	  }
	}
      }

      Eigen::VectorXd dd, cumlJK, snrval;

      if (PCAVAR) {

	// Projected coefficients
	//
	dd = (*pb)[mm]->evecJK.transpose() * (*pb)[mm]->meanJK;
      
	// Cumulative distribution
	//
	cumlJK = (*pb)[mm]->evalJK;
	for (int nn=1; nn<rank3; nn++) cumlJK[nn] += cumlJK[nn-1];
	for (int nn=0; nn<rank3; nn++) cumlJK[nn] /= cumlJK[rank3-1];
	
	// SNR vector
	//
	snrval.resize(cumlJK.size());
	
	// Compute Hall coefficients
	//
	for (int nn=0; nn<rank3; nn++) {
	  
	  double var = std::max<double>((*pb)[mm]->evalJK[nn],
					std::numeric_limits<double>::min());

	  //  Noise-to-signal ratio using the CLT estimate for
	  //  N-particle variance
	  //
	  if (subsamp)
	    var *= static_cast<double>(nbodstot) / static_cast<double>(sampT);

	  double sqr = dd[nn]*dd[nn];

	  snrval[nn] = sqrt(sqr/var);
	  
	  double b = var/sqr/nbodstot;
	  
	  (*pb)[mm]->ratio[nn] = b;

	  minSNR = std::min<double>(minSNR, 1.0/b);
	  maxSNR = std::max<double>(maxSNR, 1.0/b);
	}

#if defined HAVE_VTK
	if (vtkpca) vtkpca->Add((*pb)[mm]->meanJK,
				(*pb)[mm]->ratio, snrval,
				(*pb)[mm]->evalJK,
				(*pb)[mm]->evecJK.transpose(),
				(*pb)[mm]->covrJK,
				0, mm);
#endif
      }

      if (hout.good()) {

	double var = 0.0;

	for (int nn=0; nn<rank3; nn++) {
	  hout << setw( 4) << mm << setw(4) << nn;

	  if (PCAVAR) {
	  
	    double var = std::max<double>((*pb)[mm]->evalJK[nn],
					  std::numeric_limits<double>::min());
	    double sqr = dd[nn]*dd[nn];
	    double rat = (*pb)[mm]->ratio[nn];

	    double cof = accum_cos[mm][nn] * accum_cos[mm][nn];
	    if (mm) cof += accum_sin[mm][nn] * accum_sin[mm][nn];

	    hout << setw(fwid) << dd[nn]
		 << setw(fwid) << sqr
		 << setw(fwid) << var
		 << setw(fwid) << cumlJK[nn]
		 << setw(fwid) << snrval[nn]*snrval[nn]
		 << setw(fwid) << (*pb)[mm]->ratio[nn]
		 << setw(fwid) << 1.0/(1.0 + pow(rat, HEXP))
		 << setw(fwid) << cof
		 << setw(fwid) << (*pb)[mm]->covrJK(nn, nn)
		 << setw(fwid) << cof/(*pb)[mm]->covrJK(nn, nn)*nbodstot;
	  } else {
	    double cof = accum_cos[mm][nn] * accum_cos[mm][nn];
	    if (mm) cof += accum_sin[mm][nn] * accum_sin[mm][nn];
	    hout << setw(fwid) << sqrt(cof);
	  }
	  if (PCAEOF) hout << std::setw(fwid) << eofvec[nn];
	  hout << std::endl;
	}
	hout << std::endl;
      }
    }

#if defined HAVE_VTK
    if (vtkpca) {
      std::ostringstream sout;
      sout << hallfile << "_" << std::setfill('0') << std::setw(5) << ocount++;
      vtkpca->Write(sout.str());
    }
#endif

    // Clean storage
    //
    for (int nth=0; nth<nthrds; nth++) {

      if (PCAEOF) 
	for (auto & v : tvar[nth]) v.setZero();

      if (PCAVAR) {
	for (unsigned T=0; T<sampT; T++) {
	  numbT1[nth][T] = 0;
	  massT1[nth][T] = 0.0;
	  for (auto & v : covV[nth][T]) v.setZero();
	  for (auto & v : covM[nth][T]) v.setZero();
	}
      }
    }

    eofcount++;
  }
    
  if (VFLAG & 4)
    cerr << "Process " << setw(4) << myid << ": exiting to pca_hall" << endl;
}


void EmpCylSL::get_trimmed
(double snr,
 std::vector<Eigen::VectorXd>& ac_cos, std::vector<Eigen::VectorXd>& ac_sin,
 std::vector<Eigen::VectorXd>* rt_cos, std::vector<Eigen::VectorXd>* rt_sin,
 std::vector<Eigen::VectorXd>* sn_rat)
{
  if (PCAVAR and tk_type != None) {

    if (pb==0) return;

    ac_cos.resize(MMAX+1);
    ac_sin.resize(MMAX+1);

    if (rt_cos) {
      rt_cos->resize(MMAX+1);
      rt_sin->resize(MMAX+1);
      sn_rat->resize(MMAX+1);
    }

    // Loop through each harmonic subspace [EVEN cosines]
    //
    
    Eigen::VectorXd ddc(1, rank3), dds(1, rank3), wrk(1, rank3), sig(1, rank3);
    
    for (int mm=0; mm<=MMAX; mm++) {
      
      ac_cos[mm].resize(NORDER);
      if (rt_cos) {
	(*rt_cos)[mm].resize(NORDER);
	(*rt_sin)[mm].resize(NORDER);
	(*sn_rat)[mm].resize(NORDER);
      }
      if (mm) {
	ac_sin[mm].resize(NORDER);
      }

      sig.setZero();

      auto it = pb->find(mm);
      
      if (it != pb->end()) {
	
	auto & I = it->second;
	
	
	// Smooth coefficients
	//
	auto smth = I->ratio;

	if (tk_type == Hall) {
	  for (int i=0; i<smth.size(); i++)
	    smth[i] = 1.0/(1.0 + pow(snr*smth[i], HEXP));
	}
	if (tk_type == Truncate) {
	  for (int i=0; i<smth.size(); i++) {
	    if (1.0/smth[i]>snr) smth[i] = 1.0;
	    else                 smth[i] = 0.0;
	  }
	}


	// Cosine coefficients
	{
	  // Project to decorrelated basis
	  for (int nn=0; nn<rank3; nn++) wrk[nn] = accum_cos[mm][nn];
	  
	  ddc = I->evecJK.transpose() * wrk;

	  // Smooth
	  for (int nn=0; nn<rank3; nn++) wrk[nn] = ddc[nn] * smth[nn];

	  // Deproject coefficients
	  ddc = I->evecJK * wrk;
	  for (int nn=0; nn<rank3; nn++) ac_cos[mm][nn] = ddc[nn];

	  // Accumulate variance for modulus by error propagation
	  for (int nn=0; nn<rank3; nn++) {
	    for (int oo=0; oo<rank3; oo++) {
	      // Columns are eigenvectors
	      double val = I->evecJK(oo, nn);
	      sig[oo] += val*val*I->evalJK[nn];
	    }
	  }
	}

	// Sine coefficients
	if (mm) {
	  // Project to decorrelated basis
	  for (int nn=0; nn<rank3; nn++) wrk[nn] = accum_sin[mm][nn];
	  
	  dds = I->evecJK.transpose() * wrk;

	  // Smooth
	  for (int nn=0; nn<rank3; nn++) wrk[nn] = dds[nn] * smth[nn];

	  // Deproject coefficients
	  dds = I->evecJK * wrk;
	  for (int nn=0; nn<rank3; nn++) ac_sin[mm][nn] = dds[nn];

	} else {
	  dds.setZero();
	}

	// BEG: diagnostics
	if (rt_cos) {

	  (*rt_cos)[mm] = accum_cos[mm];
	  if (mm) (*rt_sin)[mm] = accum_sin[mm];
	  
	  for (int nn=0; nn<rank3; nn++) {
	    double val = ddc[nn]*ddc[nn] + dds[nn]*dds[nn];
	    if (sig[nn]>0.0) 
	      (*sn_rat)[mm][nn] = val/sig[nn];
	    else
	      (*sn_rat)[mm][nn] = 0.0;
	  }
	}
	// END: diagnostics
      }
      // END: covar exists
    }
    // M loop
  }
      
}

void EmpCylSL::set_trimmed(double snr, double rem)
{
  if (PCAVAR and tk_type != None) {

    if (pb==0) return;

    // Loop through each harmonic subspace [EVEN cosines]
    //
    
    Eigen::VectorXd ddc(1, rank3), dds(1, rank3), wrk(1, rank3);
    
    for (int mm=0; mm<=MMAX; mm++) {
      
      auto it = pb->find(mm);
      
      if (it != pb->end()) {
	
	auto & I = it->second;
	
	// Smooth coefficients
	//
	auto smth = I->ratio;

	if (tk_type == Hall) {
	  for (int i=0; i<smth.size(); i++)
	    smth[i] = 1.0/(1.0 + pow(snr*smth[i], HEXP));
	}
	if (tk_type == Truncate) {
	  for (int i=0; i<smth.size(); i++) {
	    if (1.0/smth[i]>snr) smth[i] = 1.0;
	    else                 smth[i] = 0.0;
	  }
	}

	// Cosine coefficients
	//
	{
	  // Project to decorrelated basis
	  for (int nn=0; nn<rank3; nn++) wrk[nn] = accum_cos[mm][nn];
	  ddc = I->evecJK.transpose() * wrk;

	  // Smooth
	  for (int nn=0; nn<rank3; nn++) wrk[nn] = ddc[nn] * smth[nn];

	  // Deproject coefficients
	  ddc = I->evecJK * wrk;
	  for (int nn=0; nn<rank3; nn++) accum_cos[mm][nn] = ddc[nn];
	}

	// Sine coefficients
	//
	if (mm) {
	  // Project to decorrelated basis
	  for (int nn=0; nn<rank3; nn++) wrk[nn] = accum_sin[mm][nn];
	  dds = I->evecJK.transpose() * wrk;

	  // Smooth
	  for (int nn=0; nn<rank3; nn++) wrk[nn] = dds[nn] * smth[nn];

	  // Deproject coefficients
	  dds = I->evecJK * wrk;
	  for (int nn=0; nn<rank3; nn++) accum_sin[mm][nn] = dds[nn];

	} else {
	  dds.setZero();
	}
      }

      if (rem>0.0) {
	double sum = 0.0, cum = 0.0;;
	for (int nn=0; nn<rank3; nn++) sum += (*pb)[mm]->evalJK[nn];
	int nf;
	for (nf=0; nf<rank3; nf++) {
	  cum += (*pb)[mm]->evalJK[nf];
	  if (1.0 - cum/sum <= rem) break;
	}
	
	// Project to decorrelated basis
	for (int nn=0; nn<rank3; nn++) wrk[nn] = accum_cos[mm][nn];
	
	dds = (*pb)[mm]->evecJK.transpose() * wrk;

	// Apply nullity
	for (int n=nf; n<rank3; n++) dds[n] = 0.0;

	// Deproject coefficients
	wrk = (*pb)[mm]->evecJK * dds;
	for (int nn=0; nn<rank3; nn++) accum_cos[mm][nn] = dds[nn];

	if (mm) {
	  
	  // Project to decorrelated basis
	  for (int nn=0; nn<rank3; nn++) wrk[nn] = accum_sin[mm][nn];
	
	  dds = (*pb)[mm]->evecJK.transpose() * wrk;

	  // Apply nullity
	  for (int n=nf; n<rank3; n++) dds[n] = 0.0;

	  // Deproject coefficients
	  wrk = (*pb)[mm]->evecJK * dds;
	  for (int nn=0; nn<rank3; nn++) accum_sin[mm][nn] = dds[nn];
	}
      }
      // END: rem
      
    }
    // M loop
  }
      
}


void EmpCylSL::accumulated_eval(double r, double z, double phi, 
				double &p0, double& p, 
				double& fr, double& fz, double &fp)
{
  if (!coefs_made_all()) {
    if (VFLAG>3)
      cerr << "Process " << myid << ": in EmpCylSL::accumlated_eval, "
	   << "calling make_coefficients()" << endl;
    make_coefficients();
  }

  fr = 0.0;
  fz = 0.0;
  fp = 0.0;
  p  = 0.0;

  double rr = sqrt(r*r + z*z);
  if (rr/ASCALE>Rtable) {
#ifdef OFF_GRID_ALERT
    std::cout << "EmpCylSL::accumulated_eval: off grid with rr=" << rr << " > "
	      << Rtable*ASCALE << std::endl;
#endif
    return;
  }

  double X = (r_to_xi(r) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;
  
  if (ix < 0) {
    ix = 0;
    if (enforce_limits) X = 0.0;
  }
  if (iy < 0) {
    iy = 0;
    if (enforce_limits) Y = 0.0;
  }
  
  if (ix >= NUMX) {
    ix = NUMX-1;
    if (enforce_limits) X = NUMX;
  }
  if (iy >= NUMY) {
    iy = NUMY-1;
    if (enforce_limits) Y = NUMY;
  }

  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;
  
  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;
  
  double ccos, ssin=0.0, fac;
  
  for (int mm=std::max<int>(0, MMIN); mm<=std::min<int>(MLIM, MMAX); mm++) {
    
    // Suppress odd M terms?
    if (EVEN_M && (mm/2)*2 != mm) continue;

    ccos = cos(phi*mm);
    ssin = sin(phi*mm);

    for (int n=std::max<int>(0, NMIN); n<std::min<int>(NLIM, rank3); n++) {
      
      fac = accum_cos[mm][n] * ccos;
      
      p += fac *
	(
	 potC[mm][n](ix  , iy  ) * c00 +
	 potC[mm][n](ix+1, iy  ) * c10 +
	 potC[mm][n](ix  , iy+1) * c01 +
	 potC[mm][n](ix+1, iy+1) * c11 
	 );
      
      fr += fac *
	(
	 rforceC[mm][n](ix  , iy  ) * c00 +
	 rforceC[mm][n](ix+1, iy  ) * c10 +
	 rforceC[mm][n](ix  , iy+1) * c01 +
	 rforceC[mm][n](ix+1, iy+1) * c11
	 );
      
      fz += fac *
	(
	 zforceC[mm][n](ix  , iy  ) * c00 +
	 zforceC[mm][n](ix+1, iy  ) * c10 +
	 zforceC[mm][n](ix  , iy+1) * c01 +
	 zforceC[mm][n](ix+1, iy+1) * c11 
	 );
      
      fac = accum_cos[mm][n] * ssin;
      
      fp += fac * mm *
	(
	 potC[mm][n](ix  , iy  ) * c00 +
	 potC[mm][n](ix+1, iy  ) * c10 +
	 potC[mm][n](ix  , iy+1) * c01 +
	 potC[mm][n](ix+1, iy+1) * c11 
	 );
      
      
      if (mm) {
	
	fac = accum_sin[mm][n] * ssin;
	
	p += fac *
	  (
	   potS[mm][n](ix  , iy  ) * c00 +
	   potS[mm][n](ix+1, iy  ) * c10 +
	   potS[mm][n](ix  , iy+1) * c01 +
	   potS[mm][n](ix+1, iy+1) * c11 
	   );
	
	fr += fac *
	  (
	   rforceS[mm][n](ix  , iy  ) * c00 +
	   rforceS[mm][n](ix+1, iy  ) * c10 +
	   rforceS[mm][n](ix  , iy+1) * c01 +
	   rforceS[mm][n](ix+1, iy+1) * c11
	   );
	
	fz += fac *
	  (
	   zforceS[mm][n](ix  , iy  ) * c00 +
	   zforceS[mm][n](ix+1, iy  ) * c10 +
	   zforceS[mm][n](ix  , iy+1) * c01 +
	   zforceS[mm][n](ix+1, iy+1) * c11 
	   );
	
	fac = -accum_sin[mm][n] * ccos;
	
	fp += fac * mm *
	  (
	   potS[mm][n](ix  , iy  ) * c00 +
	   potS[mm][n](ix+1, iy  ) * c10 +
	   potS[mm][n](ix  , iy+1) * c01 +
	   potS[mm][n](ix+1, iy+1) * c11 
	   );
	
      }
      
    }
    
    if (mm==0) p0 = p;

  }
  
}


double EmpCylSL::accumulated_dens_eval(double r, double z, double phi, 
				       double& d0)
{
  if (!DENS) return 0.0;

  if (!coefs_made_all()) {
    if (VFLAG>3) 
      cerr << "Process " << myid << ": in EmpCylSL::accumlated_dens_eval, "
	   << "calling make_coefficients()" << endl;
    make_coefficients();
  }

  double ans = 0.0;

  double rr = sqrt(r*r + z*z);

  if (rr/ASCALE > Rtable) return ans;

  double X = (r_to_xi(r) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;

  if (ix < 0) {
    ix = 0;
    if (enforce_limits) X = 0.0;
  }
  if (iy < 0) {
    iy = 0;
    if (enforce_limits) Y = 0.0;
  }
  
  if (ix >= NUMX) {
    ix = NUMX-1;
    if (enforce_limits) X = NUMX;
  }
  if (iy >= NUMY) {
    iy = NUMY-1;
    if (enforce_limits) Y = NUMY;
  }

  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;

  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;

  double ccos, ssin=0.0, fac;

  for (int mm=std::max<int>(0, MMIN); mm<=std::min<int>(MLIM, MMAX); mm++) {

    ccos = cos(phi*mm);
    ssin = sin(phi*mm);

    for (int n=std::max<int>(0, NMIN); n<std::min<int>(NLIM, rank3); n++) {

      fac = accum_cos[mm][n]*ccos;

      ans += fac *
	(
	 densC[mm][n](ix  , iy  ) * c00 +
	 densC[mm][n](ix+1, iy  ) * c10 +
	 densC[mm][n](ix  , iy+1) * c01 +
	 densC[mm][n](ix+1, iy+1) * c11 
	 );

      if (mm) {

	fac = accum_sin[mm][n]*ssin;

	ans += fac *
	  (
	   densS[mm][n](ix  , iy  ) * c00 +
	   densS[mm][n](ix+1, iy  ) * c10 +
	   densS[mm][n](ix  , iy+1) * c01 +
	   densS[mm][n](ix+1, iy+1) * c11 
	   );
      }

    }

    if (mm==0) d0 = ans;

  }

  return ans;
}


  
void EmpCylSL::get_pot(Eigen::MatrixXd& Vc, Eigen::MatrixXd& Vs,
		       double r, double z)
{
  Vc.resize(max(1,MMAX)+1, rank3);
  Vs.resize(max(1,MMAX)+1, rank3);

  if (z/ASCALE > Rtable) z =  Rtable*ASCALE;
  if (z/ASCALE <-Rtable) z = -Rtable*ASCALE;

  double X = (r_to_xi(r) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;

  if (ix < 0) {
    ix = 0;
    if (enforce_limits) X = 0.0;
  }
  if (iy < 0) {
    iy = 0;
    if (enforce_limits) Y = 0.0;
  }
  
  if (ix >= NUMX) {
    ix = NUMX-1;
    if (enforce_limits) X = NUMX;
  }
  if (iy >= NUMY) {
    iy = NUMY-1;
    if (enforce_limits) Y = NUMY;
  }

  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;

  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;

  double fac = 1.0;

  for (int mm=0; mm<=std::min<int>(MLIM, MMAX); mm++) {
    
    // Suppress odd M terms?
    if (EVEN_M && (mm/2)*2 != mm) continue;

    for (int n=0; n<rank3; n++) {

      Vc(mm, n) = fac *
	(
	 potC[mm][n](ix  , iy  ) * c00 +
	 potC[mm][n](ix+1, iy  ) * c10 +
	 potC[mm][n](ix  , iy+1) * c01 +
	 potC[mm][n](ix+1, iy+1) * c11 
	 );

      if (mm) {

	Vs(mm, n) = fac *
	  (
	   potS[mm][n](ix  , iy  ) * c00 +
	   potS[mm][n](ix+1, iy  ) * c10 +
	   potS[mm][n](ix  , iy+1) * c01 +
	   potS[mm][n](ix+1, iy+1) * c11 
	   );
      }

    }

  }

}
    
  
void EmpCylSL::get_all(int mm, int nn, 
		       double r, double z, double phi,
		       double& p, double& d, 
		       double& fr, double& fz, double &fp)
{
  if (!coefs_made_all()) {
    if (VFLAG & 4) 
      cerr << "Process " << myid << ": in EmpCylSL::gel_all, "
	   << "calling make_coefficients()" << endl;
    make_coefficients();
  }

  fr = 0.0;
  fz = 0.0;
  fp = 0.0;
  p  = 0.0;
  d  = 0.0;

  double rr = sqrt(r*r + z*z);

  if (rr/ASCALE>Rtable) {
    p = -cylmass/(rr+1.0e-16);
    fr = p*r/(rr+1.0e-16)/(rr+1.0e-16);
    fz = p*z/(rr+1.0e-16)/(rr+1.0e-16);

    return;
  }

  if (z/ASCALE > Rtable) z =  Rtable*ASCALE;
  if (z/ASCALE <-Rtable) z = -Rtable*ASCALE;

  double X = (r_to_xi(r) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;

  if (ix < 0) ix = 0;
  if (iy < 0) iy = 0;
  
  if (ix >= NUMX) ix = NUMX-1;
  if (iy >= NUMY) iy = NUMY-1;

  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;

  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;

  double ccos, ssin;

  ccos = cos(phi*mm);
  ssin = sin(phi*mm);

  p += ccos *
    (
     potC[mm][nn](ix  , iy  ) * c00 +
     potC[mm][nn](ix+1, iy  ) * c10 +
     potC[mm][nn](ix  , iy+1) * c01 +
     potC[mm][nn](ix+1, iy+1) * c11 
     );

  fr += ccos *
    (
     rforceC[mm][nn](ix  , iy  ) * c00 +
     rforceC[mm][nn](ix+1, iy  ) * c10 +
     rforceC[mm][nn](ix  , iy+1) * c01 +
     rforceC[mm][nn](ix+1, iy+1) * c11
     );
  
  fz += ccos *
    (
     zforceC[mm][nn](ix  , iy  ) * c00 +
     zforceC[mm][nn](ix+1, iy  ) * c10 +
     zforceC[mm][nn](ix  , iy+1) * c01 +
     zforceC[mm][nn](ix+1, iy+1) * c11 
     );
  
  fp += ssin * mm *
    (
     potC[mm][nn](ix  , iy  ) * c00 +
     potC[mm][nn](ix+1, iy  ) * c10 +
     potC[mm][nn](ix  , iy+1) * c01 +
     potC[mm][nn](ix+1, iy+1) * c11 
     );
  
  if (DENS)
  d += ccos *
    (
     densC[mm][nn](ix  , iy  ) * c00 +
     densC[mm][nn](ix+1, iy  ) * c10 +
     densC[mm][nn](ix  , iy+1) * c01 +
     densC[mm][nn](ix+1, iy+1) * c11 
     );


  if (mm) {
    
    p += ssin *
      (
       potS[mm][nn](ix  , iy  ) * c00 +
       potS[mm][nn](ix+1, iy  ) * c10 +
       potS[mm][nn](ix  , iy+1) * c01 +
       potS[mm][nn](ix+1, iy+1) * c11 
       );

    fr += ssin *
      (
       rforceS[mm][nn](ix  , iy  ) * c00 +
       rforceS[mm][nn](ix+1, iy  ) * c10 +
       rforceS[mm][nn](ix  , iy+1) * c01 +
       rforceS[mm][nn](ix+1, iy+1) * c11
       );

    fz += ssin *
      (
       zforceS[mm][nn](ix  , iy  ) * c00 +
       zforceS[mm][nn](ix+1, iy  ) * c10 +
       zforceS[mm][nn](ix  , iy+1) * c01 +
       zforceS[mm][nn](ix+1, iy+1) * c11 
	   );

    fp += -ccos * mm *
	  (
	   potS[mm][nn](ix  , iy  ) * c00 +
	   potS[mm][nn](ix+1, iy  ) * c10 +
	   potS[mm][nn](ix  , iy+1) * c01 +
	   potS[mm][nn](ix+1, iy+1) * c11 
	   );
      
    if (DENS)
    d += ssin *
      (
       densS[mm][nn](ix  , iy  ) * c00 +
       densS[mm][nn](ix+1, iy  ) * c10 +
       densS[mm][nn](ix  , iy+1) * c01 +
       densS[mm][nn](ix+1, iy+1) * c11 
       );

  }
  
}


void EmpCylSL::dump_coefs(ostream& out)
{
  out.setf(ios::scientific);

  for (int mm=0; mm<=MMAX; mm++) {

    out << setw(4) << mm;
    for (int j=0; j<rank3; j++)
      out << " " << setw(15) << accum_cos[mm][j];
    out << endl;

    if (mm) {
      out << setw(4) << mm;
      for (int j=0; j<rank3; j++)
	out << " " << setw(15) << accum_sin[mm][j];
      out << endl;
    }

  }
}

void EmpCylSL::set_coefs(int m1,
			 const Eigen::VectorXd& cos1, const Eigen::VectorXd& sin1, bool zero1)
{
  // Zero the coefficients
  //
  if (zero1) {
    for (int mm=0; mm<=MMAX; mm++) accum_cos[mm].setZero();
    for (int mm=1; mm<=MMAX; mm++) accum_sin[mm].setZero();

    coefs_made = vector<short>(multistep+1, true);
  }

  int nmin = std::min<int>(rank3, cos1.size());
  if (m1 <= MMAX) {
    for (int j=0; j<nmin; j++) accum_cos[m1][j] = cos1[j];
    if (m1) {
      nmin = std::min<int>(rank3, sin1.size());
      for (int j=0; j<nmin; j++) accum_sin[m1][j] = sin1[j];
    }
  }
}

void EmpCylSL::set_coefs(int m1,
			 const std::vector<double>& cos1,
			 const std::vector<double>& sin1, bool zero1)
{
  // Zero the coefficients
  //
  if (zero1) {
    for (int mm=0; mm<=MMAX; mm++) accum_cos[mm].setZero();
    for (int mm=1; mm<=MMAX; mm++) accum_sin[mm].setZero();

    coefs_made = vector<short>(multistep+1, true);
  }

  int nminC = std::min<int>(rank3, cos1.size());
  int nminS = std::min<int>(rank3, sin1.size());
  if (m1 <= MMAX) {
    for (int j=0; j<nminC; j++) accum_cos[m1][j] = cos1[j];
    if (m1) {
      for (int j=0; j<nminS; j++) accum_sin[m1][j] = sin1[j];
    }
  }
}

void EmpCylSL::get_coefs(int m1, Eigen::VectorXd& cos1, Eigen::VectorXd& sin1)
{
  if (m1 <= MMAX) {
    cos1 = accum_cos[m1];
    if (m1) sin1 = accum_sin[m1];
  }
}

void EmpCylSL::get_coefs(int m1,
			 std::vector<double>& cos1,
			 std::vector<double>& sin1)
{
  if (m1 <= MMAX) {
    cos1.resize(rank3);
    for (int j=0; j<rank3; j++) cos1[j] = accum_cos[m1][j];
    if (m1) {
      sin1.resize(rank3);
      for (int j=0; j<rank3; j++) sin1[j] = accum_sin[m1][j];
    }
  }
}


void EmpCylSL::dump_coefs_binary(ostream& out, double time)
{
  if (NewCoefs) {
    // This is a node of simple {key: value} pairs.  More general
    // content can be added as needed.
    //
    YAML::Node node;

    node["time"  ] = time;
    node["mmax"  ] = MMAX;
    node["nmax"  ] = rank3;

    // Serialize the node
    //
    YAML::Emitter y; y << node;

    // Get the size of the string
    //
    unsigned int hsize = strlen(y.c_str());
    
    // Write magic #
    //
    out.write(reinterpret_cast<const char *>(&cmagic),   sizeof(unsigned int));

    // Write YAML string size
    //
    out.write(reinterpret_cast<const char *>(&hsize),    sizeof(unsigned int));

    // Write YAML string
    //
    out.write(reinterpret_cast<const char *>(y.c_str()), hsize);

  } else {
    coefheadercyl.time = time;
    coefheadercyl.mmax = MMAX;
    coefheadercyl.nmax = rank3;

    out.write((const char *)&coefheadercyl, sizeof(CylCoefHeader));
  }

  for (int mm=0; mm<=MMAX; mm++) {

    for (int j=0; j<rank3; j++)
      out.write((const char *)&accum_cos[mm][j], sizeof(double));
    
    if (mm) {

      for (int j=0; j<rank3; j++)
	out.write((const char *)&accum_sin[mm][j], sizeof(double));
      
    }
  }
}


void EmpCylSL::dump_basis(const string& name, int step, double Rmax)
{
  static int numx = 60;
  static int numy = 60;
  
  double rmax;
  if (Rmax>0.0) rmax = Rmax;
  else rmax = 0.33*Rtable*max<double>(ASCALE, HSCALE);
  
  double r, dr = rmax/(numx-1);
  double z, dz = 2.0*rmax/(numy-1);
  double fac   = 1.0;
  int    nw    = 1 + floor(log10(0.5+NOUT));
  
  ofstream outC, outS;

  for (int mm=0; mm<=MMAX; mm++) {

    for (int n=0; n<=min<int>(NOUT, rank3-1); n++) {

      ostringstream ins;
      ins << std::setfill('0') << std::right
	  << name << ".C." << std::setw(2) << mm << "."
	  << std::setw(nw) << n << "." << step;
      
      outC.open(ins.str().c_str());
      
      if (mm) {

	ostringstream ins;
	ins << std::setfill('0') << std::right
	    << name << ".S." << std::setw(2) << mm << "."
	    << std::setw(nw) << n << "." << step;
	
	outS.open(ins.str().c_str());
      }
      
				// Ok, write data

      for (int k=0; k<numy; k++) {

	z = -rmax + dz*k;

	for (int j=0; j<numx; j++) {
	  
	  r = dr*j;

	  outC << setw(15) << r << setw(15) << z;

	  double X = (r_to_xi(r) - XMIN)/dX;
	  double Y = ( z_to_y(z) - YMIN)/dY;

	  int ix = (int)X;
	  int iy = (int)Y;

	  if (ix < 0) {
	    ix = 0;
	    if (enforce_limits) X = 0.0;
	  }
	  if (iy < 0) {
	    iy = 0;
	    if (enforce_limits) Y = 0.0;
	  }
	  
	  if (ix >= NUMX) {
	    ix = NUMX-1;
	    if (enforce_limits) X = NUMX;
	  }
	  if (iy >= NUMY) {
	    iy = NUMY-1;
	    if (enforce_limits) Y = NUMY;
	  }
	  
	  double delx0 = (double)ix + 1.0 - X;
	  double dely0 = (double)iy + 1.0 - Y;
	  double delx1 = X - (double)ix;
	  double dely1 = Y - (double)iy;

	  double c00 = delx0*dely0;
	  double c10 = delx1*dely0;
	  double c01 = delx0*dely1;
	  double c11 = delx1*dely1;

	  outC << setw(15) 
	       << fac*(
		       potC[mm][n](ix  , iy  ) * c00 +
		       potC[mm][n](ix+1, iy  ) * c10 +
		       potC[mm][n](ix  , iy+1) * c01 +
		       potC[mm][n](ix+1, iy+1) * c11 )
	       << setw(15)
	       << fac*(
		       rforceC[mm][n](ix  , iy  ) * c00 +
		       rforceC[mm][n](ix+1, iy  ) * c10 +
		       rforceC[mm][n](ix  , iy+1) * c01 +
		       rforceC[mm][n](ix+1, iy+1) * c11 )
	       << setw(15)
	       << fac*(
		       zforceC[mm][n](ix  , iy  ) * c00 +
		       zforceC[mm][n](ix+1, iy  ) * c10 +
		       zforceC[mm][n](ix  , iy+1) * c01 +
		       zforceC[mm][n](ix+1, iy+1) * c11 );
	  
	  if (DENS)
	    outC << setw(15)
		 << fac*(
			 densC[mm][n](ix  , iy  ) * c00 +
			 densC[mm][n](ix+1, iy  ) * c10 +
			 densC[mm][n](ix  , iy+1) * c01 +
			 densC[mm][n](ix+1, iy+1) * c11 );

	  outC << endl;

	  if (mm) {
      
	    outS << setw(15) << r << setw(15) << z
		 << setw(15)
		 << fac*(
			 potS[mm][n](ix  , iy  ) * c00 +
			 potS[mm][n](ix+1, iy  ) * c10 +
			 potS[mm][n](ix  , iy+1) * c01 +
			 potS[mm][n](ix+1, iy+1) * c11 )
		 << setw(15)
		 << fac*(
			 rforceS[mm][n](ix  , iy  ) * c00 +
			 rforceS[mm][n](ix+1, iy  ) * c10 +
			 rforceS[mm][n](ix  , iy+1) * c01 +
			 rforceS[mm][n](ix+1, iy+1) * c11 )
		 << setw(15)
		 << fac*(
			 zforceS[mm][n](ix  , iy  ) * c00 +
			 zforceS[mm][n](ix+1, iy  ) * c10 +
			 zforceS[mm][n](ix  , iy+1) * c01 +
			 zforceS[mm][n](ix+1, iy+1) * c11 );
	    
	    if (DENS)
	      outS << setw(15) 
		   << fac*(
			   densS[mm][n](ix  , iy  ) * c00 +
			   densS[mm][n](ix+1, iy  ) * c10 +
			   densS[mm][n](ix  , iy+1) * c01 +
			   densS[mm][n](ix+1, iy+1) * c11 );
	    outS << endl;
	  }

	}
	outC << endl;
	if (mm) outS << endl;
      }
      
				// Close output streams
      outC.close();      
      if (mm) outS.close();
    }
  }

}

void EmpCylSL::dump_images(const string& OUTFILE,
			   double XYOUT, double ZOUT, int OUTR, int OUTZ,
			   bool logscale)
{
  if (myid!=0) return;
  if (!coefs_made_all()) return;
  
  double p, d, rf, zf, pf;
  double dr, dz = 2.0*ZOUT/(OUTZ-1);
  double rmin = RMIN*ASCALE;
  
  if (logscale) 
    dr = (log(XYOUT) - log(rmin))/(OUTR-1);
  else
    dr = (XYOUT - rmin)/(OUTR-1);
  
  string Name;
  int Number     = 14;
  int Number2    = 8;
  string Types[] = {".pot", ".dens", ".fr", ".fz", ".fp",
				// Finite difference
		    ".empfr", ".empfz", ".empfp",
				// Compare with recursion
		    ".diffr", ".diffz", ".diffp",
				// Relative difference
		    ".relfr", ".relfz", ".relfp"};
  
  //============
  // Open files
  //============
  std::ofstream out[Number];
  for (int j=0; j<Number; j++) {
    Name = OUTFILE + Types[j] + ".eof_recon";
    out[j].open(Name.c_str());
    if (out[j].fail()) {
      cerr << "Couldn't open <" << Name << ">" << endl;
      break;
    }
  }
  
  double del, tmp, rP, rN, zP, zN, pP, pN, p0, d0;
  double potpr, potnr, potpz, potnz, potpp, potnp;

  
  double hr = HFAC*dr;
  double hz = HFAC*dz;
  double hp = 2.0*M_PI/64.0;
  
  if (logscale) hr = HFAC*0.5*(exp(dr) - exp(-dr));
  
  double r=0.0, z=0.0;
  
  for (int iz=0; iz<OUTZ; iz++) {
    
    z = -ZOUT + dz*iz;

    for (int ir=0; ir<OUTR; ir++) {
      
      if (logscale)
	r = rmin*exp(dr*ir);
      else
	r = rmin + dr*ir;
	  
      accumulated_eval(r, z, 0.0, p0, p, rf, zf, pf);
      d = accumulated_dens_eval(r, z, 0.0, d0);
      
      out[0]  << setw(15) << r << setw(15) << z << setw(15) << p;
      out[1]  << setw(15) << r << setw(15) << z << setw(15) << d;
      out[2]  << setw(15) << r << setw(15) << z << setw(15) << rf;
      out[3]  << setw(15) << r << setw(15) << z << setw(15) << zf;
      out[4]  << setw(15) << r << setw(15) << z << setw(15) << pf;

      
      //===================
      // Finite difference
      //===================
      
      if (logscale) {
	rP = r*(1.0 + 0.5*hr);
	rN = max<double>(1.0e-5, r*(1.0-0.5*hr));
	del = rP - rN;
      }
      else {
	rP = dr*ir+0.5*hr;
	rN = max<double>(1.0e-5, dr*ir-0.5*hr);
	del = hr;
      }
      zP = -ZOUT + dz*iz + 0.5*hz;
      zN = -ZOUT + dz*iz - 0.5*hz;
      
      pP =  0.5*hp;
      pN = -0.5*hp;

      accumulated_eval(rP, z, 0.0, tmp, potpr, tmp, tmp, tmp);
      accumulated_eval(rN, z, 0.0, tmp, potnr, tmp, tmp, tmp);
      accumulated_eval(r, zP, 0.0, tmp, potpz, tmp, tmp, tmp);
      accumulated_eval(r, zN, 0.0, tmp, potnz, tmp, tmp, tmp);
      accumulated_eval(r,  z,  pP, tmp, potpp, tmp, tmp, tmp);
      accumulated_eval(r,  z,  pN, tmp, potnp, tmp, tmp, tmp);
      
      //==================================================

      out[5] << setw(15) << r << setw(15) << z 
	     << setw(15) << -(potpr - potnr)/del;
      
      out[6] << setw(15) << r << setw(15) << z << setw(15)
	     << -(potpz - potnz)/hz;
      
      out[7] << setw(15) << r << setw(15) << z << setw(15)
	     << -(potpp - potnp)/hp;
      
      //==================================================
      

      out[8] << setw(15) << r << setw(15) << z 
	     << setw(15) << (potnr - potpr)/del - rf;
      
      out[9] << setw(15) << r << setw(15) << z << setw(15)
	     << (potnz - potpz)/hz - zf;
      
      out[10] << setw(15) << r << setw(15) << z << setw(15)
	      << (potpp - potnp)/hp - pf;
      
      //==================================================
	  
	  
      out[11] << setw(15) << r << setw(15) << z 
	      << setw(15) << ((potnr - potpr)/del - rf)/(fabs(rf)+1.0e-18);
      
      out[12] << setw(15) << r << setw(15) << z << setw(15)
	      << ((potnz - potpz)/hz - zf)/(fabs(zf)+1.0e-18);
      
      out[13] << setw(15) << r << setw(15) << z << setw(15)
	      << ((potpp - potnp)/hp - pf)/(fabs(pf)+1.0e-18);
      
      //==================================================

      for (int j=0; j<Number; j++) out[j] << endl;

    }
    
    for (int j=0; j<Number; j++) out[j] << endl;
    
  }
  
  //
  // Close current files
  //
  for (int j=0; j<Number; j++) out[j].close();

  //
  // Open files (face on)
  //
  for (int j=0; j<Number2; j++) {
    Name = OUTFILE + Types[j] + ".eof_recon_face";
    out[j].open(Name.c_str());
    if (!out[j]) {
      cerr << "Couldn't open <" << Name << ">" << endl;
      break;
    }
  }
  
  double rr, pp, xx, yy, dxy = 2.0*XYOUT/(OUTR-1);
  
  for (int iy=0; iy<OUTR; iy++) {
    yy = -XYOUT + dxy*iy;
    
    for (int ix=0; ix<OUTR; ix++) {
      xx = -XYOUT + dxy*ix;
      
      rr = sqrt(xx*xx + yy*yy);
      pp = atan2(yy, xx);
      
      accumulated_eval(rr, 0.0, pp, p0, p, rf, zf, pf);
      d = accumulated_dens_eval(rr, 0.0, pp, d0);
      
      out[0]  << setw(15) << xx << setw(15) << yy << setw(15) << p;
      out[1]  << setw(15) << xx << setw(15) << yy << setw(15) << d;
      out[2]  << setw(15) << xx << setw(15) << yy << setw(15) << rf;
      out[3]  << setw(15) << xx << setw(15) << yy << setw(15) << zf;
      out[4]  << setw(15) << xx << setw(15) << yy << setw(15) << pf;

      //===================
      // Finite difference
      //===================
      
      if (logscale) {
	rP = rr*(1.0 + 0.5*hr);
	rN = max<double>(1.0e-5, rr*(1.0-0.5*hr));
	del = rP - rN;
      }
      else {
	rP = rr+0.5*hr;
	rN = max<double>(1.0e-5, rr-0.5*hr);
	del = hr;
      }
      zP =  0.5*hz;
      zN = -0.5*hz;
      
      pP = pp + 0.5*hp;
      pN = pp - 0.5*hp;

      accumulated_eval(rP, 0.0, pp, tmp, potpr, tmp, tmp, tmp);
      accumulated_eval(rN, 0.0, pp, tmp, potnr, tmp, tmp, tmp);
      accumulated_eval(rr,  zP, pp, tmp, potpz, tmp, tmp, tmp);
      accumulated_eval(rr,  zN, pp, tmp, potnz, tmp, tmp, tmp);
      accumulated_eval(rr, 0.0, pP, tmp, potpp, tmp, tmp, tmp);
      accumulated_eval(rr, 0.0, pN, tmp, potnp, tmp, tmp, tmp);
      
      //==================================================

      out[5] << setw(15) << xx << setw(15) << yy 
	     << setw(15) << -(potpr - potnr)/del;
      
      out[6] << setw(15) << xx << setw(15) << yy << setw(15)
	     << -(potpz - potnz)/hz;
      
      out[7] << setw(15) << xx << setw(15) << yy << setw(15)
	     << -(potpp - potnp)/hp;
      
      //==================================================

      for (int j=0; j<Number2; j++) out[j] << endl;

    }
    
    for (int j=0; j<Number2; j++) out[j] << endl;
    
  }
  
  //
  // Close current files and delete
  //
  for (int j=0; j<Number2; j++) {
    try {
      out[j].close();
    }
    catch (const ofstream::failure& e) {
      std::cout << "EmpCylSL: exception closing file <"
		<< OUTFILE + Types[j] + ".eof_recon>"
		<< ": " << e.what() << std::endl;
    }
  }
}


void EmpCylSL::dump_images_basis(const string& OUTFILE,
				 double XYOUT, double ZOUT, 
				 int OUTR, int OUTZ, bool logscale,
				 int M1, int M2, int N1, int N2)
{
  if (myid!=0) return;
  
  double p, d, rf, zf, pf;
  double dr, dz = 2.0*ZOUT/(OUTZ-1);
  double rmin = RMIN*ASCALE;
  
  if (logscale) 
    dr = (log(XYOUT) - log(rmin))/(OUTR-1);
  else
    dr = (XYOUT - rmin)/(OUTR-1);
  
  int Number  = 10;
  string Types[] = {".pot", ".dens", ".fr", ".fz", ".empfr", ".empfz", 
		    ".diffr", ".diffz", ".relfr", ".relfz"};
  
  double tmp, rP, rN, zP, zN, r, z, del;
  double potpr, potnr, potpz, potnz;
  
  ofstream *out = new ofstream [Number];
  
  double hr = HFAC*dr;
  double hz = HFAC*dz;
  
  if (logscale) hr = HFAC*0.5*(exp(dr) - exp(-dr));
  
  for (int m=M1; m<=M2; m++) {
    
    for (int n=N1; n<=N2; n++) {
      
      
      //============
      // Open files
      //============
      
      for (int j=0; j<Number; j++) {
	ostringstream fname;
	fname << OUTFILE << Types[j] << '.' << m << '.' << n << ".eof_basis";
	out[j].open(fname.str().c_str());
	if (out[j].bad()) {
	  cout << "EmpCylSL::dump_images_basis: failed to open " 
	       << fname.str() << endl;
	  return;
	}
      }
      
      for (int iz=0; iz<OUTZ; iz++) {
	for (int ir=0; ir<OUTR; ir++) {
	  
	  z = -ZOUT + dz*iz;
	  if (logscale)
	    r = rmin*exp(dr*ir);
	  else
	    r = rmin + dr*ir;
	  
	  get_all(m, n, r, z, 0.0, p, d, rf, zf, pf);
	  
	  out[0] << setw(15) << r << setw(15) << z << setw(15) << p;
	  out[1] << setw(15) << r << setw(15) << z << setw(15) << d;
	  out[2] << setw(15) << r << setw(15) << z << setw(15) << rf;
	  out[3] << setw(15) << r << setw(15) << z << setw(15) << zf;
	  
	  //===================
	  // Finite difference
	  //===================
	  
	  if (logscale) {
	    rP = r*(1.0 + 0.5*hr);
	    rN = max<double>(1.0e-5, r*(1.0-0.5*hr));
	    del = rP - rN;
	  }
	  else {
	    rP = dr*ir+0.5*hr;
	    rN = max<double>(1.0e-5, dr*ir-0.5*hr);
	    del = hr;
	  }
	  zP = -ZOUT + dz*iz + 0.5*hz;
	  zN = -ZOUT + dz*iz - 0.5*hz;
	  
	  get_all(m, n, rP, z, 0.0, potpr, tmp, tmp, tmp, tmp);
	  get_all(m, n, rN, z, 0.0, potnr, tmp, tmp, tmp, tmp);
	  get_all(m, n, r, zP, 0.0, potpz, tmp, tmp, tmp, tmp);
	  get_all(m, n, r, zN, 0.0, potnz, tmp, tmp, tmp, tmp);
	  
	  out[4] << setw(15) << r << setw(15) << z 
		 << setw(15) << -(potpr - potnr)/del ;
	  
	  out[5] << setw(15) << r << setw(15) << z << setw(15)
		 << -(potpz - potnz)/hz  << setw(15) << hz;
	  
	  out[6] << setw(15) << r << setw(15) << z << setw(15)
		 << -(potpr - potnr)/del - rf;
	  
	  out[7] << setw(15) << r << setw(15) << z << setw(15)
		 << -(potpz - potnz)/hz  - zf;
	  
	  
	  out[8] << setw(15) << r << setw(15) << z << setw(15)
		 << (-(potpr - potnr)/del - rf)/(fabs(rf)+1.0e-18);
	  
	  out[9] << setw(15) << r << setw(15) << z << setw(15)
		 << (-(potpz - potnz)/hz  - zf)/(fabs(zf)+1.0e-18);
	  
	  for (int j=0; j<Number; j++) out[j] << endl;
	  
	}
	
	for (int j=0; j<Number; j++) out[j] << endl;
	
      }
      
      for (int j=0; j<Number; j++) out[j].close();
    }
    
  }
  
  delete [] out;
}

double EmpCylSL::r_to_xi(double r)
{
  if (CMAPR>0) {
    if (r<0.0) {
      ostringstream msg;
      msg << "radius=" << r << " < 0! [mapped]";
      throw GenericError(msg.str(), __FILE__, __LINE__);
    }
    return (r/ASCALE - 1.0)/(r/ASCALE + 1.0);
  } else {
    if (r<0.0)  {
      ostringstream msg;
      msg << "radius=" << r << " < 0!";
      throw GenericError(msg.str(), __FILE__, __LINE__);
    }
    return r;
  }
}
    
double EmpCylSL::xi_to_r(double xi)
{
  if (CMAPR>0) {
    if (xi<-1.0) throw GenericError("xi < -1!", __FILE__, __LINE__);
    if (xi>=1.0) throw GenericError("xi >= 1!", __FILE__, __LINE__);

    return (1.0 + xi)/(1.0 - xi) * ASCALE;
  } else {
    return xi;
  }

}

double EmpCylSL::d_xi_to_r(double xi)
{
  if (CMAPR>0) {
    if (xi<-1.0) throw GenericError("xi < -1!", __FILE__, __LINE__);
    if (xi>=1.0) throw GenericError("xi >= 1!", __FILE__, __LINE__);

    return 0.5*(1.0 - xi)*(1.0 - xi) / ASCALE;
  } else {
    return 1.0;
  }
}

#define MINEPS 1.0e-10


void EmpCylSL::legendre_R(int lmax, double x, Eigen::MatrixXd& p)
{
  /* return the associated Legendre functions up to lmax, evaluated at
     x = cos(theta) \in [-1,1].

  Using the standard forward column method from Holmes & Featherstone (2002), Section 2.1. 
  The novelty is an un-normalised associated legendre function, which remains within the limits
    of double precision to very high order.


   */
  double u,pll,norm;
  int l1,l2; // for recursion

  // start with the diagonal: l=m
  // 

  u = sqrt(1-x*x); // sin( acos(x));

  // precompute the initial values
  p(0, 0) = 1.0;
  p(1, 1) = pll = -sqrt(3)*u;

  // follow the Colombo 1981 recursion, eq. 13 of F&H2002
  if (lmax > 0) {
    for (int m=2; m<=lmax; m++) {

      // apply the recursion, eq. 13, with extra negative for Condon-Shortley convention.
      pll *= -u*sqrt((2.*m+1)/(2.*m));
      
      p(m, m) = pll;

      // check for NaNs
      if (std::isnan(p(m, m)))
	cerr << "legendre_R: p[" << m << "][" << m << "]: pllHF=" << pll << endl;
    }
  }

  // now march down the m degrees for each l order
  // (see Fig 2 of Holmes & Featherstone 2002)
  for (int l=1; l<=lmax; l++) { // start at 1 because we already have p(0,0) computed
    
    for (int m=0; m<l; m++) { // don't include m=l because we already have the diagonal p(l,l) computed

      // check that we aren't requesting an invalid location
      l1 = max(0,l - 1);
      l2 = max(0,l - 2);

      // apply the recursion, eq. 11 (and using equations 12)
      p(l,m) = sqrt((2.*l - 1)*(2*l + 1)/(l-m)/(l+m))*x*p(l1,m)
	       - sqrt((2.*l + 1)*(l+m-1)*(l-m-1)/(l-m)/(l+m)/(2*l-3))*p(l2,m);

      // version using broken-out defintions
      //p(l,m) = aml(l,m)*x*p(l1,m) - bml(l,m)*p(l2,m);
    }
  }

  
  if (std::isnan(x))
    cerr << "legendre_R: x" << endl;

  // recheck whole table for NaNs, apply overall norm to match conventions
  for(int l=0; l<=lmax; l++)
    for (int m=0; m<=l; m++) {
      if (m>0) {
	norm = sqrt(8.*M_PI);
      } else {
	norm = sqrt(4.*M_PI);
      };
      p(l,m) = p(l,m)/norm;
      
      if (std::isnan(p(l, m)))
	cerr << "legendre_R: p[" << l << "][" << m << "] lmax=" 
	     << lmax << endl;
    }
}


void EmpCylSL::dlegendre_R(int lmax, double x,
			   Eigen::MatrixXd &p, Eigen::MatrixXd &dp)
{
  /* return the associated Legendre functions up to lmax, evaluated at
     x = cos(theta) \in [-1,1].
        
         AND
  
     the derivatives, P_{lm}^(1).

  Using the standard forward column method from Holmes & Featherstone (2002), Section 2.1. 
  The novelty is an un-normalised associated legendre function, which remains within the limits
    of double precision to very high order.


   */
  int l1;
  
  // call out to fill the P_{lm}^(0) table
  legendre_R(lmax, x, p);

  if (1.0-fabs(x) < MINEPS) {
    if (x>0) x =   1.0 - MINEPS;
    else     x = -(1.0 - MINEPS);
  }

  // now we are prepared to fill out the P_{lm}^(1) table
  
  // follow whole table to fill out values
  for (int l=0; l<=lmax; l++) {
    
    for (int m=0; m<=l; m++) {
  
      // check that we aren't requesting an invalid location (computing l=m=0)
      l1 = max(0,l - 1);

      // eq. 15 of F&H
      //dp(l,m) = ( l*x*p(l,m) - fml(l,m)*p(l1,m) )/u;

      // match Martin's norm (i.e. different u definition)
      dp(l,m) = ( l*x*p(l,m) - sqrt((l*l-m*m)*(2.*l+1)/(2.*l-1))*p(l1,m) )/(x*x - 1);

      // the sectoral diagonal (l=m) follows a special relation (eq 16), but also satisfies the recursion (eq 15), so no need for anything special.

    }

  }
      
}




void EmpCylSL::sinecosine_R(int mmax, double phi,
			    Eigen::VectorXd& c, Eigen::VectorXd& s)
{
  int m;

  c[0] = 1.0;
  s[0] = 0.0;

  c[1] = cos(phi);
  s[1] = sin(phi);

  for (m=2; m<=mmax; m++) {
    c[m] = 2.0*c[1]*c[m-1] - c[m-2];
    s[m] = 2.0*c[1]*s[m-1] - s[m-2];
  }
}


void EmpCylSL::dump_eof_file(const string& eof_file, const string& output)
{
  ifstream in(eof_file.c_str());
  if (!in) {
    cerr << "---- EmpCylSL::cache_grid: error opening input file" << endl;
    return;
  }

  ofstream out(output.c_str());
  if (!in) {
    cerr << "---- EmpCylSL::cache_grid: error opening outputfile" << endl;
    return;
  }

  int mmax, numx, numy, nmax, norder, tmp, cmap;
  bool dens=false;
  double rmin, rmax, ascl, hscl;

  in.read((char *)&mmax,   sizeof(int));
  in.read((char *)&numx,   sizeof(int));
  in.read((char *)&numy,   sizeof(int));
  in.read((char *)&nmax,   sizeof(int));
  in.read((char *)&norder, sizeof(int));
  in.read((char *)&tmp,    sizeof(int)); if (tmp) dens = true;
  in.read((char *)&cmap,   sizeof(int));
  in.read((char *)&rmin,   sizeof(double));
  in.read((char *)&rmax,   sizeof(double));
  in.read((char *)&ascl,   sizeof(double));
  in.read((char *)&hscl,   sizeof(double));
  
  out << setw(70) << setfill('-') << '-' << setfill(' ') << endl;
  out << setw(20) << left << "MMAX"   << " : " << mmax << endl;
  out << setw(20) << left << "NUMX"   << " : " << numx << endl;
  out << setw(20) << left << "NUMY"   << " : " << numy << endl;
  out << setw(20) << left << "NMAX"   << " : " << nmax << endl;
  out << setw(20) << left << "NORDER" << " : " << norder << endl;
  out << setw(20) << left << "DENS"   << " : " << std::boolalpha << dens << endl;
  out << setw(20) << left << "CMAPR"  << " : " << cmap << endl;
  out << setw(20) << left << "RMIN"   << " : " << rmin << endl;
  out << setw(20) << left << "RMAX"   << " : " << rmax << endl;
  out << setw(20) << left << "ASCALE" << " : " << ascl << endl;
  out << setw(20) << left << "HSCALE" << " : " << hscl << endl;
  out << setw(70) << setfill('-') << '-' << setfill(' ') << endl;
    
  double time;
  in.read((char *)&cylmass, sizeof(double));
  in.read((char *)&time,    sizeof(double));
  out << setw(20) << left << "cylmass" << " : " << cylmass << endl;
  out << setw(20) << left << "time"    << " : " << time    << endl;
  out << setw(70) << setfill('-') << '-' << setfill(' ') << endl;

  // Read table
  //
  int nfield = 3;
  if (DENS) nfield += 1;
  
  Dynamic3dArray<double> mat(nfield, numx+1, numy+1);
  
  for (int m=0; m<=mmax; m++) {
    
    for (int v=0; v<norder; v++) {

      for (int ix=0; ix<=numx; ix++)
	for (int iy=0; iy<=numy; iy++) {
	  in.read((char *)&mat[0][ix][iy], sizeof(double));
	}
      
      for (int ix=0; ix<=numx; ix++)
	for (int iy=0; iy<=numy; iy++)
	  in.read((char *)&mat[1][ix][iy], sizeof(double));
      
      for (int ix=0; ix<=numx; ix++)
	for (int iy=0; iy<=numy; iy++)
	  in.read((char *)&mat[2][ix][iy], sizeof(double));
      
      if (DENS) {
	for (int ix=0; ix<=numx; ix++)
	  for (int iy=0; iy<=numy; iy++)
	    in.read((char *)&mat[3][ix][iy], sizeof(double));
	
      }
      
      for (int ix=0; ix<numx; ix++) {
	for (int iy=0; iy<numy; iy++) {
	  out << left << setw(4) << m << setw(4) << v 
	      << setw(4) << ix << setw(4) << iy;
	  for (int n=0; n<nfield; n++) out << setw(16) << mat[n][ix][iy]; 
	  out << endl;
	}
	
      }
      
    }
    out << setw(70) << setfill('-') << '-' << setfill(' ') << endl;

  }

  for (int m=1; m<=mmax; m++) {
    
    for (int v=0; v<norder; v++) {
      
      for (int ix=0; ix<=numx; ix++)
	for (int iy=0; iy<=numy; iy++)
	  in.read((char *)&mat[0][ix][iy], sizeof(double));
      
      for (int ix=0; ix<=numx; ix++)
	for (int iy=0; iy<=numy; iy++)
	  in.read((char *)&mat[1][ix][iy], sizeof(double));
      
      for (int ix=0; ix<=numx; ix++)
	for (int iy=0; iy<=numy; iy++)
	  in.read((char *)&mat[2][ix][iy], sizeof(double));
      
      if (DENS) {
	for (int ix=0; ix<=numx; ix++)
	  for (int iy=0; iy<=numy; iy++)
	    in.read((char *)&mat[3][ix][iy], sizeof(double));
      }

      for (int ix=0; ix<numx; ix++) {
	for (int iy=0; iy<numy; iy++) {
	  out << left << setw(4) << m << setw(4) << v 
	      << setw(4) << ix << setw(4) << iy;
	  for (int n=0; n<nfield; n++) out << setw(16) << mat[n][ix][iy]; 
	  out << endl;
	}
      }
    }
      
    out << setw(70) << setfill('-') << '-' << setfill(' ') << endl;

  }
    
}

void EmpCylSL::restrict_order(int n)
{
  for (int m=0; m<=MMAX; m++) {
    for (int k=n; k<NORDER; k++) {
      accum_cos[m][k] = 0.0;
      if (m>0) accum_sin[m][k] = 0.0;
    }
  }
}

void EmpCylSL::dump_images_basis_pca(const string& runtag,
				     double XYOUT, double ZOUT, 
				     int OUTR, int OUTZ, int M, int N, int K)
{
  if (myid!=0) return;
  if (pb == 0) return;
  
  double p, d, rf, zf, pf;

  double rmin = RMIN*ASCALE;
  double dR   = (XYOUT-rmin)/(OUTR - 1);
  double dZ   = 2.0*ZOUT/(OUTZ - 1);

  int Number  = 4;
  string Types[] = {".pot", ".dens", ".fr", ".fz"};
  
  std::vector< std::vector<double> > dataC(Number), dataS(Number);
  for (auto & v : dataC) v.resize(OUTR*OUTZ);
  if (M) for (auto & v : dataS) v.resize(OUTR*OUTZ);

  Eigen::VectorXd PP(1, NORDER), DD(1, NORDER), RF(1, NORDER), ZF(1, NORDER);
  
  std::shared_ptr<ThreeDGrid> grid;

#ifdef HAVE_VTK
  grid = std::make_shared<VtkGrid>
    (OUTR, OUTZ, 1, rmin, XYOUT, -ZOUT, ZOUT, 0, 0);
#else
  grid = std::make_shared<TableGrid>
    (OUTR, OUTZ, 1, rmin, XYOUT, -ZOUT, ZOUT, 0, 0);
#endif
  

  for (int iz=0; iz<OUTZ; iz++) {

    for (int ir=0; ir<OUTR; ir++) {
	  
      double z = -ZOUT + dZ*iz;
      double r =  rmin + dR*ir;
      double tmp;
	  
      //! Cosine space: inner produce of original basis and ev
      //! transformation

      for (int n=0; n<NORDER; n++)
	get_all(M, n, r, z, 0.0, PP[n+1], DD[n+1], RF[n+1], ZF[n+1], tmp);
      //                    ^
      //                    |
      //                    + selects COSINE only
      
      Eigen::VectorXd tp = (*pb)[M]->evecJK.row(N);

      dataC[0][ir*OUTZ + iz] = tp.adjoint() * PP;
      dataC[1][ir*OUTZ + iz] = tp.adjoint() * DD;
      dataC[2][ir*OUTZ + iz] = tp.adjoint() * RF;
      dataC[3][ir*OUTZ + iz] = tp.adjoint() * ZF;

      //! Sine space: only compute for M>0
      if (M) {
	double phi = 0.5*M_PI/M; // Selects SINE only
	for (int n=0; n<NORDER; n++)
	  get_all(M, n, r, z, phi, PP[n+1], DD[n+1], RF[n+1], ZF[n+1], tmp);

	dataS[0][ir*OUTZ + iz] = tp.adjoint() * PP;
	dataS[1][ir*OUTZ + iz] = tp.adjoint() * DD;
	dataS[2][ir*OUTZ + iz] = tp.adjoint() * RF;
	dataS[3][ir*OUTZ + iz] = tp.adjoint() * ZF;
      }
    }
  }

  for (int i=0; i<Number; i++) grid->Add(dataC[i], Types[i]+"(cos)");
  if (M) for (int i=0; i<Number; i++) grid->Add(dataS[i], Types[i]+"(sin)");
  
  std::ostringstream sout;
  sout << runtag << "_pcabasis_" << K << "_m" << M << "_n" << N;
  grid->Write(sout.str());
}


void EmpCylSL::dump_images_basis_eof(const string& runtag,
				     double XYOUT, double ZOUT, 
				     int OUTR, int OUTZ, int M, int N, int K,
				     const Eigen::VectorXd& tp)
{
  if (myid!=0) return;
  if (pb == 0) return;
  
  double p, d, rf, zf, pf;

  double rmin = RMIN*ASCALE;
  double dR   = (XYOUT-rmin)/(OUTR - 1);
  double dZ   = 2.0*ZOUT/(OUTZ - 1);

  int Number  = 4;
  string Types[] = {".pot", ".dens", ".fr", ".fz"};
  
  std::vector< std::vector<double> > dataC(Number), dataS(Number);
  for (auto & v : dataC) v.resize(OUTR*OUTZ);
  if (M) for (auto & v : dataS) v.resize(OUTR*OUTZ);

  Eigen::VectorXd PP(1, NORDER), DD(1, NORDER), RF(1, NORDER), ZF(1, NORDER);
  
  DataGrid grid(OUTR, OUTZ, 1, rmin, XYOUT, -ZOUT, ZOUT, 0, 0);

  for (int iz=0; iz<OUTZ; iz++) {

    for (int ir=0; ir<OUTR; ir++) {
	  
      double z = -ZOUT + dZ*iz;
      double r =  rmin + dR*ir;
      double tmp;
	  
      //! Cosine space: inner produce of original basis and ev
      //! transformation

      for (int n=0; n<NORDER; n++)
	get_all(M, n, r, z, 0.0, PP[n+1], DD[n+1], RF[n+1], ZF[n+1], tmp);
      //                    ^
      //                    |
      //                    + selects COSINE only
      

      dataC[0][ir*OUTZ + iz] = tp.adjoint() * PP;
      dataC[1][ir*OUTZ + iz] = tp.adjoint() * DD;
      dataC[2][ir*OUTZ + iz] = tp.adjoint() * RF;
      dataC[3][ir*OUTZ + iz] = tp.adjoint() * ZF;

      //! Sine space: only compute for M>0
      if (M) {
	double phi = 0.5*M_PI/M; // Selects SINE only
	for (int n=0; n<NORDER; n++)
	  get_all(M, n, r, z, phi, PP[n+1], DD[n+1], RF[n+1], ZF[n+1], tmp);

	dataS[0][ir*OUTZ + iz] = tp.adjoint() * PP;
	dataS[1][ir*OUTZ + iz] = tp.adjoint() * DD;
	dataS[2][ir*OUTZ + iz] = tp.adjoint() * RF;
	dataS[3][ir*OUTZ + iz] = tp.adjoint() * ZF;
      }
    }
  }

  for (int i=0; i<Number; i++) grid.Add(dataC[i], Types[i]+"(cos)");
  if (M) for (int i=0; i<Number; i++) grid.Add(dataS[i], Types[i]+"(sin)");
  
  std::ostringstream sout;
  sout << runtag << "_pcaeof_" << K << "_m" << M << "_n" << N;
  grid.Write(sout.str());
}


void EmpCylSL::compare_basis(const EmpCylSL *p)
{
  std::map<std::string, std::vector<int> > DBdif;
  std::map<std::string, std::vector<int> > DBmax;
  
  std::vector<std::string> labs =
    {"potC", "potS", "densC", "densS",
     "rforceC", "rforceS", "zforceC", "zforceS"};
  
  for (int m=0; m<=MMAX; m++) {	// Initialize values in DB
    for (auto s : labs) {
      DBdif[s].resize(MMAX+1, 0.0);
      DBmax[s].resize(MMAX+1, 0.0);
    }
  }
  
				// Compute cosine dif and maxf
  for (int m=0; m<=MMAX; m++) {

    for (int v=0; v<rank3; v++) {

      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++) {

	  double one = potC[m][v](ix, iy);
	  double two = p->potC[m][v](ix, iy);
	  
	  double cur = DBdif["potC"][m];
	  double dif = fabs(one-two);
	  DBdif["potC"][m] = dif>cur ? dif : cur;
	  
	  cur = DBmax["potC"][m];
	  dif = fabs(one);
	  DBmax["potC"][m] = dif>cur ? dif : cur;
	}
      
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++) {

	  double one = rforceC[m][v](ix, iy);
	  double two = p->rforceC[m][v](ix, iy);
	  
	  double cur = DBdif["rforceC"][m];
	  double dif = fabs(one-two);
	  DBdif["rforceC"][m] = dif>cur ? dif : cur;
	  
	  cur = DBmax["rforceC"][m];
	  dif = fabs(one);
	  DBmax["rforceC"][m] = dif>cur ? dif : cur;
	}
      
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++) {

	  double one = zforceC[m][v](ix, iy);
	  double two = p->zforceC[m][v](ix, iy);
	  
	  double cur = DBdif["zforceC"][m];
	  double dif = fabs(one-two);
	  DBdif["zforceC"][m] = dif>cur ? dif : cur;
	  
	  cur = DBmax["zforceC"][m];
	  dif = fabs(one);
	  DBmax["zforceC"][m] = dif>cur ? dif : cur;
	}
      
      if (DENS) {
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++) {

	    double one = densC[m][v](ix, iy);
	    double two = p->densC[m][v](ix, iy);
	    
	    double cur = DBdif["densC"][m];
	    double dif = fabs(one-two);
	    DBdif["densC"][m] = dif>cur ? dif : cur;
	    
	    cur = DBmax["densC"][m];
	    dif = fabs(one);
	    DBmax["densC"][m] = dif>cur ? dif : cur;
	  }
      }
      
    }
    
  }
  // END: cosine terms
  
				// Compute cosine dif and maxf
  for (int m=1; m<=MMAX; m++) {
    
    for (int v=0; v<rank3; v++) {
      
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++) {
	  
	  double one = potS[m][v](ix, iy);
	  double two = p->potS[m][v](ix, iy);
	  
	  double cur = DBdif["potS"][m];
	  double dif = fabs(one-two);
	  DBdif["potS"][m] = dif>cur ? dif : cur;
	  
	  cur = DBmax["potS"][m];
	  dif = fabs(one);
	  DBmax["potS"][m] = dif>cur ? dif : cur;
	}
      
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++) {
	  double one = rforceS[m][v](ix, iy);
	  double two = p->rforceS[m][v](ix, iy);
	  
	  double cur = DBdif["rforceS"][m];
	  double dif = fabs(one-two);
	  DBdif["rforceS"][m] = dif>cur ? dif : cur;
	  
	  cur = DBmax["rforceS"][m];
	  dif = fabs(one);
	  DBmax["rforceS"][m] = dif>cur ? dif : cur;
	}
      
      for (int ix=0; ix<=NUMX; ix++)
	for (int iy=0; iy<=NUMY; iy++) {
	  double one = zforceS[m][v](ix, iy);
	  double two = p->zforceS[m][v](ix, iy);
	  
	  double cur = DBdif["zforceS"][m];
	  double dif = fabs(one-two);
	  DBdif["zforceS"][m] = dif>cur ? dif : cur;
	  
	  cur = DBmax["zforceS"][m];
	  dif = fabs(one);
	  DBmax["zforceS"][m] = dif>cur ? dif : cur;
	}
      
      if (DENS) {
	for (int ix=0; ix<=NUMX; ix++)
	  for (int iy=0; iy<=NUMY; iy++) {
	    double one = densS[m][v](ix, iy);
	    double two = p->densS[m][v](ix, iy);
	    
	    double cur = DBdif["densS"][m];
	    double dif = fabs(one-two);
	    DBdif["densS"][m] = dif>cur ? dif : cur;
	    
	    cur = DBmax["densS"][m];
	    dif = fabs(one);
	    DBmax["densS"][m] = dif>cur ? dif : cur;
	  }
      }
      
    }
    
  }

  std::cout << "Difference values" << std::endl
	    << "-----------------" << std::endl
	    << std::setw(10) << std::right << "Table"
	    << std::setw( 5) << std::right << "m"
	    << std::setw(18) << std::right << "Dif"
	    << std::setw(18) << std::right << "Max"
	    << std::endl
	    << std::setw(10) << std::right << "------"
	    << std::setw( 5) << std::right << "--"
	    << std::setw(18) << std::right << "------"
	    << std::setw(18) << std::right << "------"
	    << std::endl;

  for (auto v : DBdif) {
    std::cout << std::setw(10) << std::right << v.first << std::endl;
    for (int m=0; m<=MMAX; m++) {
      std::cout << std::setw(10) << ""
		<< std::setw( 5) << std::right << m
		<< std::setw(18) << std::right << v.second[m]
		<< std::setw(18) << std::right << DBmax[v.first][m]
		<< std::endl;
    }
  }
  
}


// Compute non-dimensional vertical coordinate from Z
double EmpCylSL::z_to_y(double z)
{
  if (CMAPZ==1)
    return z/(fabs(z)+std::numeric_limits<double>::min())*asinh(fabs(z/HSCALE));
  else if (CMAPZ==2)
    return z/sqrt(z*z + HSCALE*HSCALE);
  else
    return z;
}

// Compute Z from non-dimensional vertical coordinate
double EmpCylSL::y_to_z(double y)
{
  if (CMAPZ==1)
    return HSCALE*sinh(y);
  else if (CMAPZ==2) {
    if (y<-1.0) throw GenericError("y < -1!", __FILE__, __LINE__);
    if (y>=1.0) throw GenericError("y >= 1!", __FILE__, __LINE__);
    return y * HSCALE/sqrt(1.0 - y*y);
  }
  else
    return y;
}

// For measure transformation in vertical coordinate
double EmpCylSL::d_y_to_z(double y)
{
  if (CMAPZ==1)
    return HSCALE*cosh(y);
  else if (CMAPZ==2) {
    if (y<-1.0) throw GenericError("y < -1!", __FILE__, __LINE__);
    if (y>=1.0) throw GenericError("y >= 1!", __FILE__, __LINE__);
    return HSCALE*pow(1.0-y*y, -1.5);
  } else
    return 1.0;
}

// Check orthogonality for basis
//
void EmpCylSL::ortho_check(std::ostream& out)
{
  if (DENS) {

    for (int mm=0; mm<=MMAX; mm++) {
      // Header
      //
      out << std::string(60, '-') << std::endl
	  << " M=" << mm << std::endl
	  << std::string(60, '-') << std::endl;

      // Normalization:
      //            +--- Gravitational energy kernel
      //            |           +--- Aximuthal
      //            |           |
      //            v           v
      double fac = -4.0*M_PI * (2.0*M_PI) * dX * dY;
      if (mm) fac *= 0.5;

      // Compute orthogonality matrix
      //
      for (int n1=0; n1<NORDER; n1++) {

	for (int n2=0; n2<NORDER; n2++) {

	  double sumC = 0.0, sumS = 0.0;

	  for (int ix=0; ix<=NUMX; ix++) {
	    double x = XMIN + dX*ix;
	    double r = xi_to_r(x);

	    for (int iy=0; iy<=NUMY; iy++) {

	      double y = YMIN + dY*iy;

	      sumC += fac * r/d_xi_to_r(x) * d_y_to_z(y) *
		potC[mm][n1](ix, iy) * densC[mm][n2](ix, iy);

	      if (mm)
		sumS += fac * r/d_xi_to_r(x) * d_y_to_z(y) *
		  potS[mm][n1](ix, iy) * densS[mm][n2](ix, iy);
	    }
	  }

	  out << std::setw(16) << sumC << std::setw(16) << sumS;
	}
	out << std::endl;
      }
    }
  } else {
    out << "EmpCylSL::ortho_check: "
	<< "can not check orthogonality without density computation"
	<< std::endl;
  }
}


void EmpCylSL::getDensSC(int mm, int n, double R, double z,
			 double& dC, double& dS)
{

  dC = 0.0;
  dS = 0.0;

  if (not DENS) return;

  if (R/ASCALE>Rtable or mm>MMAX or n>=rank3) return;

  double X = (r_to_xi(R) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;
  
  if (ix < 0) {
    ix = 0;
    if (enforce_limits) X = 0.0;
  }
  if (iy < 0) {
    iy = 0;
    if (enforce_limits) Y = 0.0;
  }
  
  if (ix >= NUMX) {
    ix = NUMX-1;
    if (enforce_limits) X = NUMX;
  }
  if (iy >= NUMY) {
    iy = NUMY-1;
    if (enforce_limits) Y = NUMY;
  }

  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;
  
  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;
  
  dC = 
    densC[mm][n](ix  , iy  ) * c00 +
    densC[mm][n](ix+1, iy  ) * c10 +
    densC[mm][n](ix  , iy+1) * c01 +
    densC[mm][n](ix+1, iy+1) * c11 ;

  if (mm)
    dS = 
      densS[mm][n](ix  , iy  ) * c00 +
      densS[mm][n](ix+1, iy  ) * c10 +
      densS[mm][n](ix  , iy+1) * c01 +
      densS[mm][n](ix+1, iy+1) * c11 ;
}


void EmpCylSL::getPotSC(int mm, int n, double R, double z,
			double& pC, double& pS)
{

  pC = 0.0;
  pS = 0.0;

  if (R/ASCALE>Rtable or mm>MMAX or n>=rank3) return;

  double X = (r_to_xi(R) - XMIN)/dX;
  double Y = (z_to_y(z)  - YMIN)/dY;

  int ix = (int)X;
  int iy = (int)Y;
  
  if (ix < 0) {
    ix = 0;
    if (enforce_limits) X = 0.0;
  }
  if (iy < 0) {
    iy = 0;
    if (enforce_limits) Y = 0.0;
  }
  
  if (ix >= NUMX) {
    ix = NUMX-1;
    if (enforce_limits) X = NUMX;
  }
  if (iy >= NUMY) {
    iy = NUMY-1;
    if (enforce_limits) Y = NUMY;
  }

  double delx0 = (double)ix + 1.0 - X;
  double dely0 = (double)iy + 1.0 - Y;
  double delx1 = X - (double)ix;
  double dely1 = Y - (double)iy;
  
  double c00 = delx0*dely0;
  double c10 = delx1*dely0;
  double c01 = delx0*dely1;
  double c11 = delx1*dely1;
  
  pC = 
    potC[mm][n](ix  , iy  ) * c00 +
    potC[mm][n](ix+1, iy  ) * c10 +
    potC[mm][n](ix  , iy+1) * c01 +
    potC[mm][n](ix+1, iy+1) * c11 ;

  if (mm)
    pS = 
      potS[mm][n](ix  , iy  ) * c00 +
      potS[mm][n](ix+1, iy  ) * c10 +
      potS[mm][n](ix  , iy+1) * c01 +
      potS[mm][n](ix+1, iy+1) * c11 ;
}
