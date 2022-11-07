#include <cstdlib>
#include <sys/time.h>
#include <sys/resource.h>
#include <limits>

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

#include <orbit.H>
#include <massmodel.H>
#include <biorth2d.H>
#include <interp.H>

#include <RespMat.H>
#include <model3d.H>
#include <sphereSL.H>
#include <TimeSeriesCoefs.H>

#include "expand.H"

#include <UserSatWake.H>

void VectorXcdSynchronize(Eigen::VectorXcd& mat, int id);
void MatrixXcdSynchronize(Eigen::MatrixXcd& mat, int id);
void ComplexSynchronize(std::complex<double>& c, int id);

inline bool isEven (int i)
{
  if ((i%2) == 0) return true;
  else return false;
}


inline bool isOdd (int i)
{
  if ((i%2) == 1) return true;
  else return false;
}

inline IntegrationType ITOSIT(int i)
{
  switch (i) {
  case ratint:
    return ratint;
  case jacoint:
    return jacoint;
  default:
    return rombint;
  }
}

double rot_matrix(int l, int m, int n, double beta);
double plgndr(int, int, double);


UserSatWake::UserSatWake(const YAML::Node& conf) : ExternalForce(conf)
{
  id = "SatelliteWake";

  LMIN 		= 0;
  LMAX 		= 4;
  MMIN 		= 0;
  MMAX 		= 4;
  lmax 		= 4;
  nmax 		= 10;
  nfreqs 	= 1200;
  HALO_TRUNC 	= -1;
  nptsE 	= 5;
  nptsK 	= 4;
  CAUCHY 	= 0;
  RATINT 	= 0;
  PTGRID 	= 200;
  NRECS 	= 512;
  DIVERGE 	= 0;
  DIVEXPON 	= 1.0;
  OLD		= 0;
  VERBOSE 	= 0;
  HALO_TYPE 	= 3;
  SITYPE 	= ITOSIT(2);
  RMODMAX 	= 20.0;
  DELTA 	= 0.0;
  OMPI 		= -0.01;
  NUMDF 	= 100;
  RA 		= 1.0e20;
  INCLINE 	= 0.0;
  PSI 		= 0.0;
  PHIP 		= 0.0;
  NUMT		= 20;
  E		= 0.0;
  Rperi 	= 0.1;
  Rsoft 	= 0.001;
  Rfac 		= 0.1;
  Mfac 		= 2.0;
  rmin 		= -1.0;
  rmax 		= -1.0;
  scale 	= 0.06667;
  numr 		= 400;
  nint 		= 12;
  Tmax 		= 4.0;
  delT 		= 0.01;
  Toffset 	= 0.0;
  satmass 	= 0.1;
  logL          = 0.0;
  XYMAX         = 1.0;
  NUMXY         = 50;
  INFILE 	= "SLGridSph.model";
  CACHEDIR 	= "./";
  UseCache 	= true;
  RespChk       = false;
  Circ          = false;

  // Constants
  //
  I 		= std::complex<double>(0.0, 1.0);
  rtol 		= 1.0e-2;

  // Default component for center
  //
  ctr_name	= "";
  cachefile 	= string(".halo_response.") + runtag;

  initialize();
							      
  if (ctr_name.size()>0) {
				// Look for the fiducial component for
				// centering
    bool found = false;
    for (auto c : comp->components) {
      if ( !ctr_name.compare(c->name) ) {
	c0 = c;
	found = true;
      break;
      }
    }

    if (!found) {
      std::ostringstream sout;
      sout << "Process " << myid << ": can't find desired component <"
	   << ctr_name << ">";
      throw GenericError(sout.str(), __FILE__, __LINE__, 35, false);
    }

  }
  else
    c0 = NULL;


  initialize_coefficients();

  userinfo();
}

UserSatWake::~UserSatWake()
{
  // Nothing
}

void UserSatWake::userinfo()
{
  if (myid) return;		// Return if node master node

  print_divider();

  cout << "** User routine SATWAKE initialized, " ;
  cout << "Filename=" << INFILE << "  DIVERGE=" << DIVERGE
       << "  DIVEXPON=" << DIVEXPON;
  if (c0) 
    cout << ", center on component <" << ctr_name << ">";
  else
    cout << ", using inertial center";
  cout << endl << setw(50) << setfill('-') << '-' << endl << setfill(' ');
  cout << "Main parameters:" << endl;
  cout << setw(9) << "" << setw(15) << "LMIN"	<< " = " << LMIN	<< endl;
  cout << setw(9) << "" << setw(15) << "LMAX"	<< " = " << LMAX	<< endl;
  cout << setw(9) << "" << setw(15) << "MMIN"	<< " = " << MMIN	<< endl;
  cout << setw(9) << "" << setw(15) << "MMAX"	<< " = " << MMAX	<< endl;
  cout << setw(9) << "" << setw(15) << "lmax"	<< " = " << lmax	<< endl;
  cout << setw(9) << "" << setw(15) << "nmax"	<< " = " << nmax	<< endl;
  cout << setw(9) << "" << setw(15) << "nfreqs"	<< " = " << nfreqs	<< endl;
  cout << setw(9) << "" << setw(15) << "HALO_TRUNC" << " = " << HALO_TRUNC << endl;
  cout << setw(9) << "" << setw(15) << "nptsE"	<< " = " << nptsE	<< endl;
  cout << setw(9) << "" << setw(15) << "nptsK"	<< " = " << nptsK	<< endl;
  cout << setw(9) << "" << setw(15) << "CAUCHY"	<< " = " << CAUCHY	<< endl;
  cout << setw(9) << "" << setw(15) << "RATINT"	<< " = " << RATINT	<< endl;
  cout << setw(9) << "" << setw(15) << "PTGRID"	<< " = " << PTGRID	<< endl;
  cout << setw(9) << "" << setw(15) << "NRECS"	<< " = " << NRECS	<< endl;
  cout << setw(9) << "" << setw(15) << "DIVERGE"<< " = " << DIVERGE	<< endl;
  cout << setw(9) << "" << setw(15) << "DIVEXPON"<< " = " << DIVEXPON	<< endl;
  cout << setw(9) << "" << setw(15) << "OLD"	<< " = " << OLD		<< endl;
  cout << setw(9) << "" << setw(15) << "VERBOSE"<< " = " << VERBOSE	<< endl;
  cout << setw(9) << "" << setw(15) << "HALO_TYPE"<< " = " << HALO_TYPE	<< endl;
  cout << setw(9) << "" << setw(15) << "SITYPE"	<< " = " << SITYPE	<< endl;
  cout << setw(9) << "" << setw(15) << "RMODMAX"<< " = " << RMODMAX	<< endl;
  cout << setw(9) << "" << setw(15) << "DELTA"	<< " = " << DELTA	<< endl;
  cout << setw(9) << "" << setw(15) << "OMPI"	<< " = " << OMPI	<< endl;
  cout << setw(9) << "" << setw(15) << "NUMDF"	<< " = " << NUMDF	<< endl;
  cout << setw(9) << "" << setw(15) << "RA"	<< " = " << RA		<< endl;
  cout << setw(9) << "" << setw(15) << "INCLINE"<< " = " << INCLINE	<< endl;
  cout << setw(9) << "" << setw(15) << "PSI"	<< " = " << PSI		<< endl;
  cout << setw(9) << "" << setw(15) << "PHIP"	<< " = " << PHIP	<< endl;
  cout << setw(9) << "" << setw(15) << "NUMT"	<< " = " << NUMT	<< endl;
  cout << setw(9) << "" << setw(15) << "E"	<< " = " << E		<< endl;
  cout << setw(9) << "" << setw(15) << "Rperi"	<< " = " << Rperi	<< endl;
  cout << setw(9) << "" << setw(15) << "Rsoft"	<< " = " << Rsoft	<< endl;
  cout << setw(9) << "" << setw(15) << "Rfac"	<< " = " << Rfac	<< endl;
  cout << setw(9) << "" << setw(15) << "Mfac"	<< " = " << Mfac	<< endl;
  cout << setw(9) << "" << setw(15) << "rmin"	<< " = " << rmin	<< endl;
  cout << setw(9) << "" << setw(15) << "rmax"	<< " = " << rmax	<< endl;
  cout << setw(9) << "" << setw(15) << "scale"	<< " = " << scale	<< endl;
  cout << setw(9) << "" << setw(15) << "numr"	<< " = " << numr	<< endl;
  cout << setw(9) << "" << setw(15) << "nint"	<< " = " << nint	<< endl;
  cout << setw(9) << "" << setw(15) << "Tmax"	<< " = " << Tmax	<< endl;
  cout << setw(9) << "" << setw(15) << "delT"	<< " = " << delT	<< endl;
  cout << setw(9) << "" << setw(15) << "Toffset"<< " = " << Toffset	<< endl;
  cout << setw(9) << "" << setw(15) << "MASS"	<< " = " << satmass	<< endl;
  cout << setw(9) << "" << setw(15) << "logL"	<< " = " << logL	<< endl;
  cout << setw(9) << "" << setw(15) << "INFILE"	<< " = " << INFILE	<< endl;
  cout << setw(9) << "" << setw(15) << "CACHEDIR"<< " = " << CACHEDIR	<< endl;
  cout << setw(9) << "" << setw(15) << "ctrname"<< " = " << ctr_name	<< endl;
  cout << setw(9) << "" << setw(15) << "UseCache"<< " = " << UseCache	<< endl;
  cout << setw(9) << "" << setw(15) << "XYMAX"	<< " = " << XYMAX	<< endl;
  cout << setw(9) << "" << setw(15) << "NUMXY"	<< " = " << NUMXY	<< endl;
  cout << setw(9) << "" << setw(15) << "RespChk"<< " = " << RespChk	<< endl;
  cout << setw(9) << "" << setw(15) << "Circ"	<< " = " << Circ	<< endl;
  cout << setw(50) << setfill('-') << '-' << endl << setfill(' ');

  I = std::complex<double>(0.0, 1.0);
  rtol = 1.0e-2;

  ctr_name = "";		// Default component for center



  cout << endl;

  print_divider();
}

const std::set<std::string>
UserSatWake::valid_keys = {
  "LMIN",
  "LMAX",
  "MMIN",
  "MMAX",
  "lmax",
  "nmax",
  "nfreqs",
  "HALO_TRUNC",
  "nptsE",
  "nptsK",
  "CAUCHY",
  "RATINT",
  "PTGRID",
  "NRECS",
  "DIVERGE",
  "DIVEXPON",
  "OLD",
  "VERBOSE",
  "HALO_TYPE",
  "SITYPE",
  "RMODMAX",
  "DELTA",
  "OMPI",
  "NUMDF",
  "RA",
  "INCLINE",
  "PSI",
  "PHIP",
  "NUMT",
  "E",
  "Rperi",
  "Rsoft",
  "Rfac",
  "Mfac",
  "rmin",
  "rmax",
  "scale",
  "numr",
  "nint",
  "Tmax",
  "delT",
  "Toffset",
  "MASS",
  "logL",
  "INFILE",
  "CACHEDIR",
  "ctrname",
  "UseCache",
  "XYMAX",
  "NUMXY",
  "RespChk",
  "Circ"
};

void UserSatWake::initialize()
{
  // Check for unmatched keys
  //
  auto unmatched = YamlCheck(conf, valid_keys);
  if (unmatched.size())
    throw YamlConfigError("UserSatWake", "parameter", unmatched,
			  __FILE__, __LINE__);

  // Assign values from YAML
  //
  try {
    if (conf["LMIN"])           LMIN               = conf["LMIN"].as<int>();
    if (conf["LMAX"])           LMAX               = conf["LMAX"].as<int>();
    if (conf["MMIN"])           MMIN               = conf["MMIN"].as<int>();
    if (conf["MMAX"])           MMAX               = conf["MMAX"].as<int>();
    if (conf["lmax"])           lmax               = conf["lmax"].as<int>();
    if (conf["nmax"])           nmax               = conf["nmax"].as<int>();
    if (conf["nfreqs"])         nfreqs             = conf["nfreqs"].as<int>();
    if (conf["HALO_TRUNC"])     HALO_TRUNC         = conf["HALO_TRUNC"].as<int>();
    if (conf["nptsE"])          nptsE              = conf["nptsE"].as<int>();
    if (conf["nptsK"])          nptsK              = conf["nptsK"].as<int>();
    if (conf["CAUCHY"])         CAUCHY             = conf["CAUCHY"].as<int>();
    if (conf["RATINT"])         RATINT             = conf["RATINT"].as<int>();
    if (conf["PTGRID"])         PTGRID             = conf["PTGRID"].as<int>();
    if (conf["NRECS"])          NRECS              = conf["NRECS"].as<int>();
    if (conf["DIVERGE"])        DIVERGE            = conf["DIVERGE"].as<int>();
    if (conf["DIVEXPON"])       DIVEXPON           = conf["DIVEXPON"].as<int>();
    if (conf["OLD"])            OLD                = conf["OLD"].as<int>();
    if (conf["VERBOSE"])        VERBOSE            = conf["VERBOSE"].as<int>();
    if (conf["HALO_TYPE"])      HALO_TYPE          = conf["HALO_TYPE"].as<int>();
    if (conf["SITYPE"])         SITYPE             = ITOSIT(conf["SITYPE"].as<int>());
    if (conf["RMODMAX"])        RMODMAX            = conf["RMODMAX"].as<double>();
    if (conf["DELTA"])          DELTA              = conf["DELTA"].as<double>();
    if (conf["OMPI"])           OMPI               = conf["OMPI"].as<double>();
    if (conf["NUMDF"])          NUMDF              = conf["NUMDF"].as<int>();
    if (conf["RA"])             RA                 = conf["RA"].as<double>();
    if (conf["INCLINE"])        INCLINE            = conf["INCLINE"].as<double>();
    if (conf["PSI"])            PSI                = conf["PSI"].as<double>();
    if (conf["PHIP"])           PHIP               = conf["PHIP"].as<double>();
    if (conf["NUMT"])           NUMT               = conf["NUMT"].as<int>();
    if (conf["E"])              E                  = conf["E"].as<double>();
    if (conf["Rperi"])          Rperi              = conf["Rperi"].as<double>();
    if (conf["Rsoft"])          Rsoft              = conf["Rsoft"].as<double>();
    if (conf["Rfac"])           Rfac               = conf["Rfac"].as<double>();
    if (conf["Mfac"])           Mfac               = conf["Mfac"].as<double>();
    if (conf["rmin"])           rmin               = conf["rmin"].as<double>();
    if (conf["rmax"])           rmax               = conf["rmax"].as<double>();
    if (conf["scale"])          scale              = conf["scale"].as<double>();
    if (conf["numr"])           numr               = conf["numr"].as<int>();
    if (conf["nint"])           nint               = conf["nint"].as<int>();
    if (conf["Tmax"])           Tmax               = conf["Tmax"].as<double>();
    if (conf["delT"])           delT               = conf["delT"].as<double>();
    if (conf["Toffset"])        Toffset            = conf["Toffset"].as<double>();
    if (conf["MASS"])           satmass            = conf["MASS"].as<double>();
    if (conf["logL"])           logL               = conf["logL"].as<double>();
    if (conf["INFILE"])         INFILE             = conf["INFILE"].as<string>();
    if (conf["CACHEDIR"])       CACHEDIR           = conf["CACHEDIR"].as<string>();
    if (conf["ctrname"])        ctr_name           = conf["ctrname"].as<string>();
    if (conf["UseCache"])       UseCache           = conf["UseCache"].as<int>();
    if (conf["XYMAX"])          XYMAX              = conf["XYMAX"].as<double>();
    if (conf["NUMXY"])          NUMXY              = conf["NUMXY"].as<double>();
    if (conf["RespChk"])        RespChk            = conf["RespChk"].as<int>();
    if (conf["Circ"])           Circ               = conf["Circ"].as<int>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in UserSatWake: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }
}  

void UserSatWake::initialize_coefficients()
{
  //==========================
  // Begin normal exectution
  //==========================
  
  const int id  = 0x12;
  std::complex<double>       omp_halo;
  
  //=================
  // Sanity checks
  //=================
  
  if (MMIN<0 && myid==0) cerr << "Warning: Mmin should be >=0 . . . fixing\n";
  MMIN = abs(MMIN);
  if (MMIN>LMIN) {
    if (myid==0) cout << "Setting LMIN=MMIN=" << MMIN << endl;
    LMIN = MMIN;
  }
  
  
  // ===================================================================
  // Initilize HALO model
  // ===================================================================
  
  // SphericalModelTable::linear = 1;
  auto m = std::make_shared<SphericalModelTable>(INFILE, DIVERGE, DIVEXPON);
  m->setup_df(NUMDF, RA);
  halo_model = m;
  Model3dNames[0] = INFILE;	// Assign filename to ID string
  
  // ===================================================================
  // Initilize biorthogonal functions
  // ===================================================================
  
  switch (HALO_TYPE) {
  case bessel:
    u = std::make_shared<BSSphere>(RMODMAX, nmax, LMAX);
    break;
  case clutton_brock:
    u = std::make_shared<CBSphere>();
    break;
  case hernquist:
    u = std::make_shared<HQSphere>();
    break;
  case sturm:
    if (rmin<0.0) 
      rmin = halo_model->get_min_radius() * 1.2;
    if (rmax<0.0) 
      rmax = halo_model->get_max_radius() * 0.99;

    SphereSL::mpi = 1;
    u = std::make_shared<SphereSL>(LMAX, nmax, numr, rmin, rmax, scale, m);
    break;
  default:
    {
      std::ostringstream sout;
      sout << "Illegal spherical biorthongal series: " << HALO_TYPE;
      throw GenericError(sout.str(), __FILE__, __LINE__);
    }
  }
  
  
  // ===================================================================
  // Setup orbit
  // ===================================================================
  
  if (Circ)
    E = 0.5*halo_model->get_mass(Rperi)/Rperi + halo_model->get_pot(Rperi);

  double MaxOm = sqrt(halo_model->get_mass(Rfac)/(Rfac*Rfac*Rfac));
  if (myid==0) cout << "Omega(" << Rfac << ")=" << MaxOm  << endl;

  MaxOm *= Mfac*LMAX;
  
  vector< vector<RespMat> > total;

  vector<std::complex<double>> freqs(nfreqs);
  double dOmega = 2.0*MaxOm/(nfreqs-1);
  for (int k=0; k<nfreqs; k++)
    freqs[k] = std::complex<double>(-MaxOm + dOmega*k, OMPI);

  if (myid==0) {
    cout << "Max(Omega)=" << MaxOm  << endl
	 << "    NFreqs=" << nfreqs << endl;
  }

  Times = vector<double>(NUMT+1);
  for (int nt=0; nt<=NUMT; nt++) Times[nt] = -Tmax + 2.0*Tmax*nt/NUMT;
  
  string coutfile("");
  if (myid==0) coutfile = outdir + runtag;

  TimeSeriesCoefs Coefs(E, Rperi, Rsoft, delT, Tmax, halo_model, satmass, logL, 
			coutfile);
  
  //
  // This is a map of maps . . . a nice way to make sparse matrix
  //
  std::map< int, std::map<int, std::vector<Eigen::MatrixXcd> > > coefs;
  Eigen::VectorXcd tcoefs;
  int icnt;
  
  // ===================================================================
  // Get combined response matrix for disk and halo for each L
  // ===================================================================
  
  ofstream to_save;
  ifstream from_save;
  unsigned short reading = (UseCache ? 1 : 0);

  if (myid==0 && UseCache) {
    
    from_save.open(string(outdir + cachefile).c_str());

    if (!from_save) {
      cerr << "Couldn't open <" << cachefile <<
	"> to read cached data" << endl;
      cerr << ". . .  will attempt to write a cache file." << endl;
      
      reading = 0;

    } else {
      
      // Read ID string and check version #
      int tid;
      from_save.read((char *)&tid, sizeof(int));
      if (tid != id) {
	throw GenericError("Incompatible save file!",
			   __FILE__, __LINE__, -34, false);
      } 
    
      char tbuf[255];
      from_save.read(tbuf, 255);
      cout << tbuf << endl;

      int nfreq1, lmin1, lmax1, mmin1, mmax1;
      from_save.read((char *)&nfreq1, sizeof(int));
      from_save.read((char *)&lmin1,  sizeof(int));
      from_save.read((char *)&lmax1,  sizeof(int));
      from_save.read((char *)&mmin1,  sizeof(int));
      from_save.read((char *)&mmax1,  sizeof(int));

      if (nfreqs != nfreq1) reading = 0;
      if (LMIN   != lmin1 ) reading = 0;
      if (LMAX   != lmax1 ) reading = 0;
      if (MMIN   != mmin1 ) reading = 0;
      if (MMAX   != mmax1 ) reading = 0;

      if (reading) {
	std::complex<double> tmp;
	for (int i=0; i<nfreqs; i++) {
	  from_save.read((char *)&tmp, sizeof(std::complex<double>));
	  if (fabs(tmp - freqs[i])>1.0e-10) reading = 0;
	}
      }

      if (reading==0) {
	std::ostringstream sout;
	sout << "Cache file <" << outdir << cachefile
	     << "> is incompatible with current input parameters!!";
	throw GenericError(sout.str(), __FILE__, __LINE__);
      }
    }
  }
    
  MPI_Bcast(&reading, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);
  
  if (myid && reading) {

    from_save.open(string(outdir + cachefile).c_str());

    // Read ID string and check version #
    int tid;
    from_save.read((char *)&tid, sizeof(int));
    if (tid != id) {
      std::ostringstream sout;
      sout << "Process " << myid << ": error reading save file!";
      throw GenericError(sout.str(), __FILE__, __LINE__);
    } 

    char tbuf[255];
    from_save.read(tbuf, 255);

    int nfreq1, lmin1, lmax1, mmin1, mmax1;
    from_save.read((char *)&nfreq1, sizeof(int));
    from_save.read((char *)&lmin1,  sizeof(int));
    from_save.read((char *)&lmax1,  sizeof(int));
    from_save.read((char *)&mmin1,  sizeof(int));
    from_save.read((char *)&mmax1,  sizeof(int));

    std::complex<double> tmp;
    for (int i=0; i<nfreqs; i++) {
      from_save.read((char *)&tmp, sizeof(std::complex<double>));
    }
  }
    
  if (!reading && myid==0 && UseCache) {
      
    to_save.open(string(outdir + cachefile).c_str());
    if (!to_save) {
      std::ostringstream sout;
      sout << "Couldn't open <" << cachefile << "> to write cached data";
      throw GenericError(sout.str(), __FILE__, __LINE__);
    }
    
    // Write ID string and version #
    time_t tp = time(NULL);
    string stime = string("UserSatWake") + ": " + ctime(&tp);
    cout << "ID string: " << stime<< endl;
    
    to_save.write((const char *)&id, sizeof(int));
    char tbuf[255]; tbuf[254] = 0;
    strncpy(tbuf, stime.c_str(), 254);
    to_save.write(tbuf, 255);
    
    to_save.write((const char *)&nfreqs, sizeof(int));
    to_save.write((const char *)&LMIN,   sizeof(int));
    to_save.write((const char *)&LMAX,   sizeof(int));
    to_save.write((const char *)&MMIN,   sizeof(int));
    to_save.write((const char *)&MMAX,   sizeof(int));

    for (int i=0; i<nfreqs; i++) {
      to_save.write((const char *)&freqs[i], sizeof(std::complex<double>));
    }
  }


  if (myid==0) cout << "Computing Laplace coefficients . . ." << endl;
  
  //
  // Round-robin work queue
  //
  icnt = 0;
  for (int L=LMIN; L<=LMAX; L++) {

    for (int L2=-L; L2<=L; L2+=2) {

      int id = icnt++ % numprocs;
      if (id == myid) {
	Coefs.coefs(L, L2, nmax, nint, u, freqs, Times, coefs[L][L2]);
      }
    }
  }


  //
  // Share with all processes
  //
  icnt = 0;
  for (int L=LMIN; L<=LMAX; L++) {

    if (myid==0) cout << "L=" << L << endl;

    for (int L2=-L; L2<=L; L2+=2) {
      int id = icnt++ % numprocs;
      unsigned sz = coefs[L][L2].size();
      MPI_Bcast(&sz, 1, MPI_UNSIGNED, id, MPI_COMM_WORLD);

      if (id==myid) {
	for (unsigned j=0; j<sz; j++)
	  MatrixXcdSynchronize(coefs[L][L2][j], id);
      } else {
	Eigen::MatrixXcd tmp;
	for (unsigned j=0; j<sz; j++) {
	  MatrixXcdSynchronize(tmp, id);
	  coefs[L][L2].push_back(tmp);
	}
      }
      
      if (myid==0) cout << "    L2=" << L2 << " [rows, cols]=[" 
			<< coefs[L][L2][0].rows()  << ", " 
			<< coefs[L][L2][0].cols() << "]" << endl;
    }
  }

  if (myid==0) cout << "done" << endl;
  
  // Factor for point mass expansion for
  // orbit in x-y plane
  
  INCLINE *= M_PI/180.0;
  PSI     *= M_PI/180.0;
  PHIP    *= M_PI/180.0;
  
  std::complex<double> I(0.0, 1.0);
  
  for (int L=LMIN; L<=LMAX; L++) {
    int mmin;
    if (isOdd(L)) {
      if (isOdd(MMIN)) mmin = MMIN;
      else             mmin = MMIN+1;
    } else {
      if (isOdd(MMIN)) mmin = MMIN+1;
      else             mmin = MMIN;
    }
    for (int M=mmin; M<=L; M+=2) {
      Lhalo.push_back(L);
      Mhalo.push_back(M);
    }
  }

  Nhalo = Lhalo.size();
  if (myid==0) {
    cout << "Harmonics:" << endl << "---------" << endl
	 << setw(3) << left << "L"
	 << setw(3) << left << "M" << endl
	 << setw(3) << left << "-"
	 << setw(3) << left << "-" << endl;
    for (int i=0; i<Nhalo; i++)
      cout << setw(3) << left << Lhalo[i]
	   << setw(3) << left << Mhalo[i]
	   << endl;
  }
  
  icnt = 0;
  for (int nf=0; nf<nfreqs; nf++) {
    
    if (myid==0) cout << "Halo frequency: " << freqs[nf] << endl;
    
    vector<RespMat> mab_halo(Nhalo);

    if (reading) {
      
      for (int ihalo=0; ihalo<Nhalo; ihalo++)
	mab_halo[ihalo] = RespMat(from_save, halo_model, u);
      
      total.push_back(mab_halo);
      
    } else {
      
      //
      // Get response matrices for each halo l-order
      //
      for (int ihalo=0; ihalo<Nhalo; ihalo++) {
	
	mab_halo[ihalo] = RespMat(freqs[nf], Lhalo[ihalo], Mhalo[ihalo], 
				  lmax, nmax, nptsK, nptsE,
				  halo_model, u, OLD, 0, 0,
				  VERBOSE, SITYPE, CAUCHY, RATINT);
	
	mab_halo[ihalo].set_params(DELTA, PTGRID);
	
	int id = icnt++ % numprocs;

	if (id == myid) {
	  if (CAUCHY)
	    mab_halo[ihalo].make_matrix_pv();
	  else 
	    mab_halo[ihalo].make_matrix();
	}
	
      }
      
      total.push_back(mab_halo);
    }

  }

  if (!reading) {

    icnt = 0;

    for (int nf=0; nf<nfreqs; nf++) {

      for (int ihalo=0; ihalo<Nhalo; ihalo++) {

	int id = icnt++ % numprocs;
	total[nf][ihalo].MPISynchronize(id);

	if (UseCache && myid==0)
	  total[nf][ihalo].write_out(to_save);

      }

    }

  }

  
  // ===================================================================
  // For paranoia's sake
  // ===================================================================

  if (UseCache) {
    if (reading)
      from_save.close();
    else
      if (myid==0) to_save.close();
  }  

  // ===================================================================
  // Compute the time series of response vectors
  // ===================================================================
  
  std::ofstream tlog;
  if (myid==0) tlog.open((outdir + runtag + ".satwake.timelog").c_str());
  
  Eigen::VectorXcd tmp(nmax);
  tmp.setZero();

  //
  // All the time slices for each (L,M) pair
  //

  rcoefs.resize(NUMT+1);
  for (int nt=0; nt<=NUMT; nt++) {
    for (int ihalo=0; ihalo<Nhalo; ihalo++)
      rcoefs[nt].push_back(tmp);
  }
  
  //
  // For interpolation in compute_coefficients()
  //

  curcoefs = vector<Eigen::VectorXcd>(Nhalo);
  for (int ihalo=0; ihalo<Nhalo; ihalo++)
    curcoefs.push_back(tmp);


  icnt = 0;

  for (int nt=0; nt<=NUMT; nt++) {
    
    if (myid==0) tlog << setw(8) << nt << setw(16) << Times[nt] << endl;
    
    //
    // Get satellite perturbation coefficients
    //
    for (int ihalo=0; ihalo<Nhalo; ihalo++) {
      
      int id = icnt++ % numprocs;
      
      if (id != myid) continue;

      int L = Lhalo[ihalo];
      int M = Mhalo[ihalo];
      double dup = M ? 2.0 : 1.0;
	
      for (int L2=-L; L2<=L; L2+=2) {

	//
	// Frequency loop for current time
	//
	for (int nf=0; nf<nfreqs; nf++) {
	  
	  // 
	  // Truncate satellite coefficients
	  //
	  tcoefs = coefs[L][L2][nt];
	  if (HALO_TRUNC>0) {
	    for (int n=HALO_TRUNC+1; n<=nmax; n++) tcoefs[n] = 0.0;
	  }
	  
	  //
	  // Factor for norm in density component of biorthogonal pair
	  //         |
	  //         |          Account for pos and neg M values 
	  //         |          |
	  //         |          |     Factor from spherical harmonic
	  //         |          |     |
	  //         |          |     |
	  //         V          V     V
	  tcoefs *= 4.0*M_PI * dup * sqrt( (0.5*L + 0.25)/M_PI *
            exp(lgamma(1.0+L-abs(L2)) - lgamma(1.0+L+abs(L2))) ) *
            plgndr(L, abs(L2), 0.0) * 
            rot_matrix(L, M, L2, INCLINE) *
	    exp(I*static_cast<double>(L2)*PSI) *
	    exp(I*static_cast<double>(M)*PHIP) ;
	  //  ^
	  //  |
	  //  |
	  //  Rotate component to spherical harmonic in new orientation
	  //
	  
	  
	  //
	  // Leading sign for spherical harmonic with m>0 (following Edmonds
	  // to be consistent with definition of rotation matrices)
	  
	  if (L2>0 && isOdd(L2)) tcoefs *= -1.0;

	  tcoefs *= dOmega/2.0*M_PI;
	
	  rcoefs[nt][ihalo] 
	    += total[nf][ihalo].get_response(tcoefs, RespMat::self);
	  
	}
      }

      
      for (int n=1; n<=nmax; n++)
	cout << setw(4) << nt << setw(4) << ihalo << rcoefs[nt][ihalo][n]
	     << endl;
    }
  }

  
  //
  // Share with all processes
  //
  icnt = 0;

  for (int nt=0; nt<=NUMT; nt++) {
    
    for (int ihalo=0; ihalo<Nhalo; ihalo++) {
      int id = icnt++ % numprocs;
      VectorXcdSynchronize(rcoefs[nt][ihalo], id);
    }

  }

  if (RespChk) {
    check_response();
				// For debugging
    ostringstream ostr;
    ostr << outdir << runtag << "_coefs." << myid;
    ofstream outc(ostr.str().c_str());
    if (outc) {
      for (int nt=0; nt<=NUMT; nt++) {
	for (int ihalo=0; ihalo<Nhalo; ihalo++) {
	  for (int n=1; n<=nmax; n++) 
	    outc << setw(18) << Times[nt]
		 << setw( 5) << Lhalo[ihalo]
		 << setw( 5) << Mhalo[ihalo]
		 << setw( 5) << n
		 << setw(18) << rcoefs[nt][ihalo][n].real()
		 << setw(18) << rcoefs[nt][ihalo][n].imag()
		 << setw(18) << fabs(rcoefs[nt][ihalo][n])
		 << endl;
	  outc << endl;
	}
      }
    }
    
    if (myid==0) {
      for (int ihalo=0; ihalo<Nhalo; ihalo++) {
	ostringstream ostr;
	ostr << outdir << runtag 
	     << "." << Lhalo[ihalo]
	     << "." << Mhalo[ihalo]
	     << ".coefs";
	ofstream out(ostr.str().c_str());
	for (int nt=0; nt<=NUMT; nt++) {
	  for (int n=1; n<=nmax; n++) 
	    out << setw(18) << Times[nt]
		<< setw( 5) << n
		<< setw(18) << rcoefs[nt][ihalo][n].real()
		<< setw(18) << rcoefs[nt][ihalo][n].imag()
		<< setw(18) << fabs(rcoefs[nt][ihalo][n])
		<< endl;
	  out << endl;
	}
      }
    }
  }
}
    

void UserSatWake::determine_acceleration_and_potential(void)
{
#if HAVE_LIBCUDA==1		// Cuda compatibility
  getParticlesCuda(cC);
#endif

  compute_coefficients();
  exp_thread_fork(false);
}


void * UserSatWake::determine_acceleration_and_potential_thread(void * arg) 
{
  int id = *((int*)arg);
  unsigned nbodies = cC->Number();
  int nbeg = nbodies*id/nthrds;
  int nend = nbodies*(id+1)/nthrds;

  double pos[3], rr, r, x, y, z, phi, costh, pot0, pot1, dpot;
  double potr, potl, pott, potp, R2, fac, Ylm, dYlm;
  Eigen::VectorXd f, g;
  
  PartMapItr it = cC->Particles().begin();

  for (int q=0   ; q<nbeg; q++) it++;
  for (int q=nbeg; q<nend; q++) {

    unsigned long indx = (it++)->first;
				// If we are multistepping, compute accel 
				// only at or below this level

    if (multistep && (cC->Part(indx)->level < mlevel)) continue;

    rr = 0.0;
    for (int k=0; k<3; k++) {
      pos[k] = cC->Pos(indx, k);	// Inertial by default
      if (c0) pos[k] -= c0->center[k];
      rr += pos[k]*pos[k];
    }
    r = sqrt(rr);

    potr = potl = pott = potp = 0.0;

    if (r <= halo_model->get_max_radius()) {
      
      x = pos[0];
      y = pos[1];
      z = pos[2];

      phi = atan2(y, x);
      costh = z/(r+DSMALL);
      R2 = x*x + y*y;
      
      for (unsigned ihalo=0; ihalo<Lhalo.size(); ihalo++) {

	int L = Lhalo[ihalo];
	int M = Mhalo[ihalo];

	fac = sqrt( (0.5*L + 0.25)/M_PI * 
		    exp(lgamma(1.0+L-M) - lgamma(1.0+L+M)) ) * satmass;
	Ylm = fac * plgndr(L, M, costh);
	dYlm = fac * dplgndr(L, M, costh);

	double dM = M;
	f = ( exp(I*phi*dM) * curcoefs[ihalo] ).real();
	g = ( I*dM*exp(I*phi*dM) * curcoefs[ihalo] ).real();
	  
	pot0 = u->get_potl(r, L, f);
	pot1 = u->get_potl(r, L, g);
	dpot = (u->get_potl(r*(1.0+rtol), L, f) - 
		u->get_potl(r*(1.0-rtol), L, f) ) / ( 2.0*rtol*r );
	
	potl += Ylm  * pot0;
	potr += Ylm  * dpot;
	pott += dYlm * pot0;
	potp += Ylm  * pot1;
      }
      
      cC->AddAcc(indx, 0, -(potr*x/r - pott*x*z/(r*r*r)) );
      cC->AddAcc(indx, 1, -(potr*y/r - pott*y*z/(r*r*r)) );
      cC->AddAcc(indx, 2, -(potr*z/r + pott*R2 /(r*r*r))   );
      
      if (R2 > DSMALL) {
	cC->AddAcc(indx, 0,  potp*y/R2 );
	cC->AddAcc(indx, 1, -potp*x/R2 );
      }
      if (use_external)
	cC->AddPotExt(indx, potl);
      else
	cC->AddPot(indx, potl);
    }
  }

  return (NULL);
}


void UserSatWake::compute_coefficients()
{
  int indx;
  double a, b, curT = tnow - Toffset;

  //
  // Find time index
  //
  if (curT<Times.front()) {	// No perturbation if time is out of bounds
    indx = 0;
    a = 0.0;
    b = 0.0;
  } else if (curT>=Times.back()) {
    indx = Times.size()-2;
    a = 0.0;
    b = 0.0;
  } else {
				// Index bounds enforcement
    indx = max<int>(0, min<int>(Vlocate(curT, Times), Times.size()-2));
    a = (Times[indx+1] - curT)/(Times[indx+1] - Times[indx]);
    b = (curT - Times[indx])/(Times[indx+1] - Times[indx]);
  }

  for (int ihalo=0; ihalo<Nhalo; ihalo++)
    curcoefs[ihalo] = a*rcoefs[indx][ihalo] + b*rcoefs[indx+1][ihalo];

}


void UserSatWake::check_response()
{
  int icnt = 0;
  for (int nt=0; nt<=NUMT; nt++) {
    
    int id = icnt++ % numprocs;
    if (id != myid) continue;
    
    cout << "Process " << myid << ": printing T=" << Times[nt] << endl;

    ostringstream suffix;
    suffix << nt;

    string name1 = outdir + "dens." + runtag + ".abs." + suffix.str();
    string name2 = outdir + "dens." + runtag + ".rel." + suffix.str();
    string name3 = outdir + "potl." + runtag + ".abs." + suffix.str();
    string name4 = outdir + "potl." + runtag + ".rel." + suffix.str();
    //
    // Print out wakes
    //
    gnuplot_out(rcoefs[nt],
		RespMat::self, RespMat::density, false, name1);

    gnuplot_out(rcoefs[nt],
		RespMat::self, RespMat::density, true,  name2);

    gnuplot_out(rcoefs[nt],
		RespMat::self, RespMat::potential, false, name3);

    gnuplot_out(rcoefs[nt],
		RespMat::self, RespMat::potential, true,  name4);
  }

}

void UserSatWake::gnuplot_out(vector<Eigen::VectorXcd>& coefs, 
			      RespMat::gravity grav, RespMat::response type,
			      bool relative, string& file)
{
  double x, dx = 2.0*XYMAX/(NUMXY-1);
  double y, dy = 2.0*XYMAX/(NUMXY-1);
  double r, phi, ans, back;
  constexpr std::complex<double> I(0.0, 1.0);
  Eigen::VectorXd f;
  
  std::ofstream out(file.c_str());
  if (!out) {
    std::cerr << "Could not open <" << file << "> for output" << std::endl;
    return;
  }

  for (int i=0; i<NUMXY; i++) {
    x = -XYMAX + dx*i;

    for (int j=0; j<NUMXY; j++) {
      y = -XYMAX + dy*j;
      
      r = sqrt(x*x + y*y);
      
      ans = 0.0;

      if (r <= halo_model->get_max_radius()) {
	phi = atan2(y, x);
      
	for (unsigned ihalo=0; ihalo<Lhalo.size(); ihalo++) {

	  int L = Lhalo[ihalo];
	  int M = Mhalo[ihalo];
	  if (L<M) {
	    cerr << "L=" << L << " but M=" << M << endl;
	  }

	  double Ylm = satmass*
	    sqrt( (0.5*L + 0.25)/M_PI * 
		  exp(lgamma(1.0+L-M) - lgamma(1.0+L+M)) ) * plgndr(L, M, 0.0);

	  double dM = M;
	  f = ( exp(I*phi*dM) * coefs[ihalo] ).real();
	  
	  switch (type) {
	  case density:
	    ans += Ylm*u->get_dens(r, L, f)/(4.0*M_PI);
	    break;
	  case potential:
	    ans += Ylm*u->get_potl(r, L, f);
	    break;
	  }
	}
      
	if (relative) {
	  switch (type) {
	  case density:
	    if (relative) {
	      back = halo_model->get_density(r);
	      if (fabs(back)>0.0) ans /= back;
	    }
	    break;
	  case potential:
	    if (relative) {
	      back = halo_model->get_pot(r);
	      if (fabs(back)>0.0) ans /= back;
	    }
	    break;
	  }
	}
      }

      out << setw(16) << x << setw(16) << y << setw(16) << ans << endl;
    }
    out << endl;
  }
}



extern "C" {
  ExternalForce *satWake(const YAML::Node& conf)
  {
    return new UserSatWake(conf);
  }
}

class proxysatwake { 
public:
  proxysatwake()
  {
    factory["usersatwake"] = satWake;
  }
};

static proxysatwake p;
