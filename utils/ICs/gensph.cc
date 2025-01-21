/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Generate a phase space realization from a distribution function
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/20/91, ported to EXP tree 9/28/21
 *
 ***************************************************************************/

#include <filesystem>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <memory>

// EXP classes
//
#include <libvars.H>
#include <localmpi.H>
#include <numerical.H>
#include <massmodel.H>
#include <model3d.H>
#include <GenPoly.H>
#include <isothermal.H>
#include <hernquist_model.H>
#include <EllipForce.H>
#include <localmpi.H>
#include <fpetrap.h>
#include <cxxopts.H>
#include <EXPini.H>

// Reference to n-body globals
//
using namespace __EXP__;

#include <config_exp.h>

#ifdef HAVE_FFTW
#include <complex>
#include <fftw3.h>
#endif

#include <euler.H>

// Global variables

int
main(int argc, char **argv)
{
  int HMODEL, N, NUMDF, NUMR, NUMJ, NUME, NUMG, NREPORT, SEED, ITMAX;
  int NUMMODEL, RNUM, DIVERGE, DIVERGE2, LINEAR, NUMINT, NI, ND, Nangle;
  int Nrepl=1, Nfib=1;
  double DIVERGE_RFAC, DIVERGE_RFAC2, NN, MM, RA, RMODMIN, RMOD, EPS;
  double X0, Y0, Z0, U0, V0, W0, TOLE;
  double Emin0, Emax0, Kmin0, Kmax0, RBAR, MBAR, BRATIO, CRATIO, SMOOTH;
  bool LOGR, ELIMIT, VERBOSE, GRIDPOT, MODELS, EBAR, zeropos, zerovel;
  bool VTEST;
  std::string INFILE, MMFILE, OUTFILE, OUTPS, config;

  const double goldenRatio = 1.618033988749895;

#ifdef DEBUG
  set_fpu_handler();
#endif

#ifdef FPETRAP
  fpeinit(0);
#endif

  // MPI initialization
  //
  local_init_mpi(argc, argv);

  // Option parsing
  //
  cxxopts::Options options("gensph", "Generate single-mass or multi-mass spherical ICs");

  options.add_options()
    ("h,help", "Print this help message")
    ("v,verbose", "Print additional diagnostic information")
    ("T,template", "Write template options file with current and all default values")
    ("c,config", "Parameter configuration file",
     cxxopts::value<string>(config))
    ("i,INFILE", "Mass model file",
     cxxopts::value<string>(INFILE)->default_value("mass.table"))
    ("n,MMFILE", "Number model file",
     cxxopts::value<string>(MMFILE))
    ("o,PSFILE", "Phase-space output file",
     cxxopts::value<string>(OUTPS)->default_value("new.bods"))
    ("p,prefix", "Diagnostic output file prefix",
     cxxopts::value<string>(OUTFILE)->default_value("gensph"))
    ("zeropos", "Set the origin at the center of mass",
     cxxopts::value<bool>(zeropos)->default_value("false"))
    ("zerovel", "Set the total momentum of the realization to zero",
     cxxopts::value<bool>(zerovel)->default_value("true"))
    ("HMODEL", "Halo type (0=file)",
     cxxopts::value<int>(HMODEL)->default_value("0"))
    ("N,bodies", "Number of bodies",
     cxxopts::value<int>(N)->default_value("1000000"))
    ("NI", "Number of integer attributesfor PSP file",
     cxxopts::value<int>(NI)->default_value("0"))
    ("ND", "Number of double attributes for PSP file",
     cxxopts::value<int>(ND)->default_value("0"))
    ("NUMDF", "Number of points in energy grid for Eddington inversion",
     cxxopts::value<int>(NUMDF)->default_value("10000"))
    ("NUMR", "Number of points in radial grid for halo model",
     cxxopts::value<int>(NUMR)->default_value("400"))
    ("NUMJ", "Number of points in kappa grid for Eddington inversion",
     cxxopts::value<int>(NUMJ)->default_value("400"))
    ("NUME", "Number of energy points in PS generation table",
     cxxopts::value<int>(NUME)->default_value("400"))
    ("NUMG", "Number of pointsin mass table for model realization",
     cxxopts::value<int>(NUMG)->default_value("800"))
    ("NREPORT", "Report after generating NREPORT points",
     cxxopts::value<int>(NREPORT)->default_value("1000"))
    ("SEED", "Initial seed for random number generator",
     cxxopts::value<int>(SEED)->default_value("11"))
    ("ITMAX", "Maximum number of interations for acceptance-rejection method",
     cxxopts::value<int>(ITMAX)->default_value("100000"))
    ("NUMMODEL", "Number of points for GeneralizedPolytrope",
     cxxopts::value<int>(NUMMODEL)->default_value("500"))
    ("RNUM", "Number of radial points for interally computed mass model",
     cxxopts::value<int>(RNUM)->default_value("10000"))
    ("DIVERGE", "Use inner cusp extrapolation on real mass-density model",
     cxxopts::value<int>(DIVERGE)->default_value("0"))
    ("DIVERGE_RFAC", "Inner cusp slope",
     cxxopts::value<double>(DIVERGE_RFAC)->default_value("1.5"))
    ("DIVERGE2", "Use inner cusp extrapolation on (pseudo) number-density model",
     cxxopts::value<int>(DIVERGE2)->default_value("0"))
    ("DIVERGE_RFAC2", "Inner cusp slope",
     cxxopts::value<double>(DIVERGE_RFAC2)->default_value("1.5"))
    ("LOGR", "Use logarithmic mapping for internal radial grid",
     cxxopts::value<bool>(LOGR)->default_value("false"))
    ("LINEAR", "Use linear interpolation for SphericalModelTable",
     cxxopts::value<int>(LINEAR)->default_value("1"))
    ("NN", "First polytropic index (energy)",
     cxxopts::value<double>(NN)->default_value("2.5"))
    ("MM", "Second polytropic index (ang. mom.)",
     cxxopts::value<double>(MM)->default_value("0.5"))
    ("RA", "Anisotropy index",
     cxxopts::value<double>(RA)->default_value("1.0e8"))
    ("RMODMIN", "Inner radius for Istothermal and Hernquist model",
     cxxopts::value<double>(RMODMIN)->default_value("1.0e-2"))
    ("RMOD", "Outer radius for Isothermal and Hernquist model",
     cxxopts::value<double>(RMOD)->default_value("100.0"))
    ("EPS", "step size for computing polytrope",
     cxxopts::value<double>(EPS)->default_value("1.0e-5"))
    ("X0", "Phase space offset",
     cxxopts::value<double>(X0)->default_value("0.0"))
    ("Y0", "Phase space offset",
     cxxopts::value<double>(Y0)->default_value("0.0"))
    ("Z0", "Phase space offset",
     cxxopts::value<double>(Z0)->default_value("0.0"))
    ("U0", "Phase space offset",
     cxxopts::value<double>(U0)->default_value("0.0"))
    ("V0", "Phase space offset",
     cxxopts::value<double>(V0)->default_value("0.0"))
    ("W0", "Phase space offset",
     cxxopts::value<double>(W0)->default_value("0.0"))
    ("TOLE", "Point generation fractional energy offset for Eddington grid",
     cxxopts::value<double>(TOLE)->default_value("1.0e-4"))
    ("ELIMIT", "Limit particle selection with energy and kappa bounds",
     cxxopts::value<bool>(ELIMIT)->default_value("false"))
    ("VTEST", "Test gen_velocity() generation",
     cxxopts::value<bool>(VTEST)->default_value("false"))
    ("Nangle", "Number of angular points on a regular spherical grid (0 skips this algorithm).  Nangle>0 takes precidence over other algorithms.",
     cxxopts::value<int>(Nangle)->default_value("0"))
    ("Nrepl", "Number of replicates in orbital plane (1 skips the Sellwood 1997 algorithm)",
     cxxopts::value<int>(Nrepl)->default_value("1"))
    ("Nfib", "Number of points on the sphere for each orbit. Replicate the orbits by tiling the angular momentum direction on the sphere using a Fibonnaci sequence.  Default value of 1 implies one plane per orbit.",
     cxxopts::value<int>(Nfib)->default_value("1"))
    ("Emin0", "Minimum energy (if ELIMIT=true)",
     cxxopts::value<double>(Emin0)->default_value("-3.0"))
    ("Emax0", "Maximum energy (if ELIMIT=true)",
     cxxopts::value<double>(Emax0)->default_value("-1.0"))
    ("Kmin0", "Minimum kappa (if ELIMIT=true)",
     cxxopts::value<double>(Kmin0)->default_value("0.0"))
    ("Kmax0", "Maximum kappa (if ELIMIT=true)",
     cxxopts::value<double>(Kmax0)->default_value("1.0"))
    ("RBAR", "Semi-major axis for bar ellipsoid",
     cxxopts::value<double>(RBAR)->default_value("0.067"))
    ("MBAR", "Mass of bar ellipsoid",
     cxxopts::value<double>(MBAR)->default_value("0.00103739"))
    ("BRATIO", "axis ratio b/a",
     cxxopts::value<double>(BRATIO)->default_value("0.2"))
    ("CRATIO", "axis ratio c/b",
     cxxopts::value<double>(CRATIO)->default_value("0.05"))
    ("SMOOTH", "Edge smoothing for bar density profile",
     cxxopts::value<double>(SMOOTH)->default_value("0.0"))
    ("g,gridpot", "compute the number-density potential by gridded quadrature",
     cxxopts::value<bool>(GRIDPOT)->default_value("false"))
    ("d,diagout", "print model computation diagnostics for multimass and bar",
     cxxopts::value<bool>(MODELS)->default_value("false"))
    ("b,bar",     "add an ellipsoidal bar model",
     cxxopts::value<bool>(EBAR)->default_value("false"))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid == 0) {
      std::cout << std::string(80, '%') << std::endl
		<< options.help() << std::endl << std::endl;
    }

    MPI_Finalize();
    return 0;
  }

  // Write template config file in INI style and exit
  //
  if (vm.count("template")) {
    // Write YAML template file
    //
    if (myid==0) SaveConfig(vm, options, "template.yaml");

    MPI_Finalize();
    return 0;
  }

  // Read parameters fron the YAML config file
  //
  if (vm.count("config")) {
    try {
      vm = LoadConfig(options, config);
    } catch (cxxopts::OptionException& e) {
      if (myid==0) std::cout << "Option error in configuration file: "
			     << e.what() << std::endl;
      MPI_Finalize();
      return 0;
    }
  }

  if (vm.count("verbose")) VERBOSE = true;
  else                     VERBOSE = false;

  // Use Fibonnaci for space-filling grid points
  //
  if (Nangle>0) {
    Nfib = Nrepl = 1;
  }
  else {
    Nfib  = std::max<int>(1, Nfib );
    Nrepl = std::max<int>(1, Nrepl);
  }

  // Prepare output streams and create new files
  //
  std::ostringstream sout;
  sout << OUTPS << "." << myid;
  std::ofstream out(sout.str());
  int bad = 0;
  if (!out) {
    std::cerr << "[" << myid << "] Couldn't open <" << sout.str()
	      << "> for output" << std::endl;
    bad = 1;
  }
  MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (bad) {
    MPI_Finalize();
    exit(-1);
  }

  out.precision(11);
  out.setf(ios::scientific, ios::floatfield);

  // Begin integration
  //
  AxiSymModPtr hmodel;
  SphModTblPtr htmodel, htmodel2, mmmodel;
  SphModMultPtr multi;

  SphericalModelTable::even = 0;

  AxiSymModel::numr     = NUMR;
  AxiSymModel::numj     = NUMJ;
  AxiSymModel::gen_N    = NUMG;
  AxiSymModel::gen_E    = NUME;
  AxiSymModel::gen_tolE = TOLE;

  if (HMODEL>=0) {

    // Halo model
    //
    switch (HMODEL) {
    case file:
      SphericalModelTable::even = 0;
      SphericalModelTable::linear = LINEAR;
      htmodel = std::make_shared<SphericalModelTable>(INFILE, DIVERGE, DIVERGE_RFAC);

      htmodel->setup_df(NUMDF, RA);

      if (VERBOSE and myid==0) {
	htmodel->print_df("gensph.df");
      }

      hmodel = htmodel;
      break;
    case isothermal:
      hmodel = std::make_shared<IsothermalSphere>(1.0, RMODMIN, RMOD);
      break;
    case hernquist_model:
      hmodel = std::make_shared<HernquistSphere>(1.0, RMODMIN, RMOD);
      break;
    case gen_polytrope:
      hmodel = std::make_shared<GeneralizedPolytrope>(NUMMODEL, NN, MM, EPS);
      break;
    default:
      if (myid==0)
	std::cerr << "No such HALO model type: " << HMODEL << std::endl;
      exit(-2);
    }

  }


  double rmin = hmodel->get_min_radius();
  double rmax = hmodel->get_max_radius();
  double r, dr = (log(rmax) - log(rmin))/(RNUM-1);

  vector<double> r2(RNUM), d2(RNUM), m2(RNUM), p2(RNUM), t2(RNUM), MS(RNUM);

  if (LOGR)
    dr = (log(rmax) - log(rmin))/(RNUM-1);
  else
    dr = (rmax - rmin)/(RNUM-1);


  for (int i=0; i<RNUM; i++) {
    if (LOGR)
      r2[i] = rmin*exp(dr*i);
    else
      r2[i] = rmin + dr*i;
  }

  //--------------------
  // Make bar if desired
  //--------------------

  if (EBAR) {

    if (myid==0) cout << "-----------" << std::endl;
    auto ellip = std::make_shared<EllipForce>
      (RBAR, RBAR*BRATIO, RBAR*BRATIO*CRATIO, MBAR, NUMINT, 200);

#ifdef HAVE_FFTW
    if (SMOOTH>0.0) {

      std::vector<double> a(NUMR), b(NUMR), c(NUMR);
      std::vector<std::complex<double>> A(NUMR), B(NUMR), C(NUMR);
      fftw_plan pa, pb, pinv;

      pa   = fftw_plan_dft_r2c_1d(NUMR, a.data(),
				  reinterpret_cast<fftw_complex*>(A.data()),
				  FFTW_ESTIMATE);

      pb   = fftw_plan_dft_r2c_1d(NUMR, b.data(),
				  reinterpret_cast<fftw_complex*>(B.data()),
				  FFTW_ESTIMATE);

      pinv = fftw_plan_dft_c2r_1d(NUMR,
				  reinterpret_cast<fftw_complex*>(C.data()),
				  c.data(), FFTW_ESTIMATE);

      double xmin=rmin, xmax=rmax;

      double x, z, dx = (30*SMOOTH + xmax - xmin)/NUMR;
      double y, dk = (30.0*SMOOTH + xmax - xmin)/NUMR;
      double scale = dk / NUMR;
      int indx;

      std::ofstream tin, tout;

      if (myid==0) {
	tin.open("ebar_fft.input");
	tout.open("ebar_fft.output");
      }

      for (int i=0; i<NUMR; i++) {

	x = xmin + dx*i;
	indx = (NUMR/2+i) % NUMR;
	z = dx*(indx-NUMR/2);

	if (x<xmax)
	  a[indx] = ellip->getMass(x);
	else
	  a[indx] = 0.0;

	b[indx] = exp( -z*z/(2.0*SMOOTH*SMOOTH)) / sqrt(2.0*M_PI*SMOOTH*SMOOTH);

	if (myid==0)
	  tin << std::setw(10) << indx
	      << std::setw(20) << x
	      << std::setw(20) << a[indx]
	      << std::setw(20) << z
	      << std::setw(20) << b[indx]
	      << std::endl;
      }

      fftw_execute(pa);
      fftw_execute(pb);

      for (int i=0; i<NUMR; i++) {
	C[i] = A[i] * B[i] * scale;
      }

      // inverse transform to get c, the convolution of a and b;
      //            this has the side effect of overwriting C
      fftw_execute(pinv);

      // ...

      double mbar=-1.0;
      std::vector<double> rr(NUMR), mm(NUMR);
      double fac;

      for (int i=0; i<NUMR; i++) {

	rr[i] = x = xmin + dx*i;
	y = ellip->getMass(x);
	fac = 0.5*(1.0 + erf((x-0.1*RBAR)/(0.025*RBAR)));
	if (rr[i]>RBAR+30.0*SMOOTH && mbar<0.0) mbar = MS[i-1];
	if (mbar>0.0)
	  mm[i] = mbar;
	else
	  mm[i] = (1.0 - fac)*y + fac*c[i-1];

	if (myid==0)
	  tout << std::setw(20) << x
	       << std::setw(20) << y
	       << std::setw(20) << mm[i]
	       << std::endl;
      }

      fftw_destroy_plan(pa);
      fftw_destroy_plan(pb);
      fftw_destroy_plan(pinv);

      Linear1d mass1(rr, mm);

      for (int i=0; i<RNUM; i++) MS[i] = mass1.eval(r2[i]);

    } else {
      for (int i=0; i<RNUM; i++) MS[i] = ellip->getMass(r2[i]);
    }
#else
    for (int i=0; i<RNUM; i++) MS[i] = ellip->getMass(r2[i]);
#endif

    //
    // Make new model including dark halo and bar
    //

    for (int i=0; i<RNUM; i++) {
      d2[i] = hmodel->get_density(r2[i]);
      m2[i] = hmodel->get_mass(r2[i]) + MS[i];
    }


    double dm1, dm2;
    p2[0] = 0.0;
    t2[0] = 0.0;
    dm2 = odd2(r2[0], r2, m2);
    for (int i=1; i<RNUM; i++) {
      dm1 = dm2;
      dm2 = drv2(r2[i], r2, m2);
      t2[i] += 0.5*(dm1/r2[i-1] + dm2/r2[i])*(r2[i] - r2[i-1]) + t2[i-1];
    }

    for (int i=0; i<RNUM; i++) {
      if (r2[i]>0.0)
	p2[i] = - m2[i]/r2[i] - (t2[RNUM-1] - t2[i]);
      else
	p2[i] = - (t2[RNUM-1] - t2[i]);
    }


    htmodel2 = std::make_shared<SphericalModelTable>
      (RNUM, r2.data(), d2.data(), m2.data(), p2.data());

    htmodel2->setup_df(NUMDF, RA);

    hmodel = htmodel2;

    if (MODELS and myid==0) {
      string mname = OUTFILE + ".ellip";
      ofstream outmod(mname.c_str());
      if (outmod) {
	outmod.setf(ios::scientific);
	outmod.precision(11);
	outmod << "# Internal model with ellip" << std::endl;
	outmod << RNUM << std::endl;
	for (int i=0; i<RNUM; i++) {
	  if (LOGR)
	    r = rmin*exp(dr*i);
	  else
	    r = rmin + dr*i;

	  outmod << std::setw(20) << r2[i]
		 << std::setw(20) << d2[i]
		 << std::setw(20) << m2[i]
		 << std::setw(20) << p2[i]
		 << std::endl;
	}
      } else {
	cerr << "Error opening <" << mname << std::endl;
      }
    }

  }

  int nfib=0, nangle=0, nshell=0;
  if (Nangle) {
    nangle  = Nangle;
    nfib    = floor(sqrt(double(N)/(nangle*nangle)));
    nshell  = nfib*nangle;

    // Shells in radius and velocity
    //
    N = nshell*nshell;
    if (myid==0)
    std::cout << std::setw(60) << std::setfill('-') << '-'
	      << std::setfill(' ') << std::endl
	      << "Fibonacci: N=" << N << " Nangle=" << Nangle
	      << " Nmag=" << nfib << " Nshell=" << nshell<< std::endl;
  }

  // For replication algorithm; Nrepl=1 is the standard algorithm.
  //
  Nrepl *= Nfib;
  int nplane = N/Nrepl;
  N = Nrepl * nplane;
  NREPORT = std::max<int>(1, NREPORT/Nrepl);

  if (Nrepl > 1 and myid==0) {
    std::cout << std::setw(60) << std::setfill('-') << '-'
	      << std::setfill(' ') << std::endl
	      << "Replication: N=" << N << " Nrepl=" << Nrepl/Nfib
	      << " Nfib=" << Nfib << " Norbits=" << nplane << std::endl;
  }

  double mass = hmodel->get_mass(hmodel->get_max_radius())/N;
  AxiSymModPtr rmodel = hmodel;

  if (vm.count("MMFILE")) {

    SphericalModelTable::even = 0;
    SphericalModelTable::linear = LINEAR;

    // Generate "fake" profile
    //
    mmmodel = std::make_shared<SphericalModelTable>
      (MMFILE, DIVERGE2, DIVERGE_RFAC2);

    mass = mmmodel->get_mass(mmmodel->get_max_radius())/N;

    // Packs fake density and mass model with target (real) potential
    // and reinitializes the model
    //
    std::vector<double> r2(RNUM), d2(RNUM), m2(RNUM), p2(RNUM);

    double rlast = rmin;
    vector<double> mt2(RNUM, 0.0), tt2(RNUM, 0.0);

    if (GRIDPOT) {
      m2[0] = 0.0;
      for (int i=1; i<RNUM; i++) {
	if (LOGR)
	  r2[i] = r = rmin*exp(dr*i);
	else
	  r2[i] = r = rmin + dr*i;

	mt2[i] = mt2[i-1] +
	  0.5*(hmodel->get_density(rlast)*rlast*rlast +
	       hmodel->get_density(r2[i])*r2[i]*r2[i] )
	  * (r2[i] - rlast) * 4.0*M_PI;

	m2[i] = m2[i-1] +
	  0.5*(mmmodel->get_density(rlast)*rlast*rlast +
	       mmmodel->get_density(r2[i])*r2[i]*r2[i] )
	  * (r2[i] - rlast) * 4.0*M_PI;

	tt2[i] = tt2[i-1] +
	  0.5*(hmodel->get_density(rlast)*rlast +
	       hmodel->get_density(r2[i])*r2[i] )
	  * (r2[i] - rlast) * 4.0*M_PI;

	rlast = r2[i];
      }

      if (VERBOSE and myid==0)
	cout << "-----------" << std::endl
	     << "Mass (computed)=" << mt2[RNUM-1] << std::endl;
    }

    for (int i=0; i<RNUM; i++) {
      if (LOGR)
	r2[i] = r = rmin*exp(dr*i);
      else
	r2[i] = r = rmin + dr*i;

      d2[i] = mmmodel -> get_density(r);

      if (GRIDPOT) {
	if (r2[i]<=0.0)
	  p2[i] = 0.0;
	else
	  p2[i] = -mt2[i]/r2[i] - (tt2[RNUM-1] - tt2[i]);
      } else {
	m2[i] = mmmodel -> get_mass(r);
	p2[i] = hmodel  -> get_pot(r);
      }

    }

    mmmodel = std::make_shared<SphericalModelTable>
      (RNUM, r2.data(), d2.data(), m2.data(), p2.data(),
       DIVERGE2, DIVERGE_RFAC2);

    mmmodel->setup_df(NUMDF, RA);

    if (VERBOSE and myid==0) {
      mmmodel->print_df("gensph.dfmulti");
    }

    if (MODELS and myid==0) {
      string mname = OUTFILE + ".multi";
      ofstream outmod(mname.c_str());
      if (outmod) {
	outmod.setf(ios::scientific);
	outmod.precision(11);
	outmod << "# Internal multimass model" << std::endl;
	outmod << RNUM << std::endl;
	for (int i=0; i<RNUM; i++) {
	  if (LOGR)
	    r = rmin*exp(dr*i);
	  else
	    r = rmin + dr*i;

	  outmod << std::setw(20) << r
		 << std::setw(20) << mmmodel->get_density(r)
		 << std::setw(20) << mmmodel->get_mass(r)
		 << std::setw(20) << mmmodel->get_pot(r)
		 << std::endl;
	}
      } else {
	cerr << "Error opening <" << mname << std::endl;
      }
    }

    // Generate the multimass model
    //
    if (EBAR)
      multi = std::make_shared<SphericalModelMulti>(htmodel2, mmmodel);
    else
      multi = std::make_shared<SphericalModelMulti>(htmodel,  mmmodel);

    if (vm.count("allow")) multi->allowNegativeMass();

    rmodel = multi;
  }

  random_gen.seed(SEED*numprocs+myid);
  rmodel->set_itmax(ITMAX);

  if (MODELS and myid==0) {
    std::ofstream outmod(OUTFILE);
    if (outmod) {
      outmod.setf(ios::scientific);
      outmod.precision(11);
      outmod << "# Internal model" << std::endl;
      outmod << RNUM << std::endl;
      for (int i=0; i<RNUM; i++) {
	if (LOGR)
	  r = rmin*exp(dr*i);
	else
	  r = rmin + dr*i;

	outmod << std::setw(20) << r
	       << std::setw(20) << rmodel->get_density(r)
	       << std::setw(20) << rmodel->get_mass(r)
	       << std::setw(20) << rmodel->get_pot(r)
	       << std::endl;
      }
    } else {
      std::cerr << "Error opening <" << OUTFILE << std::endl;
    }

  }

  int ierr;
  Eigen::VectorXd ps(7), ps0(7);

  ps0[0] = 0.0;
  ps0[1] = X0;
  ps0[2] = Y0;
  ps0[3] = Z0;
  ps0[4] = U0;
  ps0[5] = V0;
  ps0[6] = W0;

  if (myid==0) {
    out << std::setw(12) << N
	<< std::setw( 6) << NI << std::setw( 6) << ND  << std::endl;
  }

  // Diagnostic variables
  //
  int count=0, negms=0;
  double TT=0.0, WW=0.0, VC=0.0;

  if (myid==0)
    std::cout << std::setw(60) << std::setfill('-') << '-' << std::endl
	      << std::setfill('-')
	      << "Body count:" << std::endl;

  int npernode = nplane/numprocs;
  int beg = myid*npernode;
  int end = beg + npernode;

  if (myid==numprocs-1) end = nplane;

  std::vector<Eigen::VectorXd> PS;
  Eigen::VectorXd zz = Eigen::VectorXd::Zero(7);

  std::vector<double> rmass;
  if (Nangle) {
    rmass.resize(nfib);

    // Set up mass grid
    double rmin = hmodel->get_min_radius();
    double rmax = hmodel->get_max_radius();
    double mmas = hmodel->get_mass(rmax);
    double dmas = mmas/nfib;
    for (int i=0; i<nfib; i++) {
      double mas = dmas*(0.5+i);
      auto loc = [&](double r) -> double
      { return (mas - hmodel->get_mass(r)); };
      rmass[i] = zbrent(loc, rmin, rmax, 1.0e-8);
    }
  }

  for (int n=beg; n<end; n++) {

    do {
      if (Nangle) {
	// Get space and velocity indices
	int nspace = n/nshell;
	int nveloc = n - nspace*nshell;

	// Break space index into radius and angle
	int nr     = nspace/nangle;
	int rm     = nspace - nr*nangle;

	// Sanity check
	assert((nspace>=0 && nspace<nshell));
	assert((nveloc>=0 && nveloc<nshell));
	assert((rm>=0 && rm<nangle));

	// Get spatial coordinates
	double r    = rmass[nr];
	double phi  = 2.0*M_PI * rm / goldenRatio;
	double cost = 1.0 - 2.0 * rm / nangle;
	
	double PHI   = 2.0*M_PI * nr / goldenRatio;
	double COST  = 1.0 - 2.0 * nr / nangle;
	double THETA = acos(COST) + 0.5*M_PI;

	// Get velocity coorindates
	std::tie(ps, ierr) = rmodel->gen_point_iso(r, acos(cost), phi,
						   nveloc, nfib, nangle,
						   PHI, THETA, PHI);

	for (int i=0; i<3; i++) {
	  if (std::isnan(ps[i+3])) {
	    std::cout << "Vel NaN with r=" << r << " cost=" << cost << " phi=" << phi << std::endl;
	  }
	}
      }
      else if (ELIMIT)
	std::tie(ps, ierr) = rmodel->gen_point(Emin0, Emax0, Kmin0, Kmax0);
      else if (VTEST) {
	std::tie(ps, ierr) = rmodel->gen_point();
	if (ierr==0) {
	  ierr = rmodel->gen_velocity(&ps[1], &ps[4]);
	  if (ierr) {
	    std::cout << "gen_velocity failed: "
		      << ps[0] << " "
		      << ps[1] << " "
		      << ps[2] << "\n";
	  }
	}
      }
      else
	std::tie(ps, ierr) = rmodel->gen_point();
      if (ierr) count++;
    } while (ierr);

    if (ps[0] <= 0.0) negms++;

    Eigen::Vector3d pos(ps[1], ps[2], ps[3]), pos1;
    Eigen::Vector3d vel(ps[4], ps[5], ps[6]), vel1;
    Eigen::Matrix3d proj, Iprj, rot = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d projM, IprjM;

    double dq = 2.0*M_PI * Nfib / Nrepl, mfac = ps[0];

    // Compute orbital plane basis
    //
    if (Nrepl>1) {
      auto L = pos.cross(vel);	// The angular momentum vector
      if (pos.norm() < 1.0e-10 or L.norm() < 1.0e-10) {
	proj = Iprj = Eigen::Matrix3d::Identity();
      } else {
	auto X = pos/pos.norm(); // Radial unit vector, X.
	auto Y = L.cross(X);	 // Choose a right-hand3d coordinate
	Y = Y/Y.norm();		 // system perpendicular to X and L.
	auto Z = L/L.norm();     // Unit vector in the ang. mom. direction.
	
	// Azimuth of vertical plane containing L
	double Phi   = atan2(Z(1), Z(0)) + 0.5*M_PI;

	// Tilt of plane w.r.t to polar axis
	double Theta = acos(Z(2));

	// Location of X in orbital plane w.r.t line of nodes
	Eigen::Vector3d lon(cos(Phi), sin(Phi), 0.0);
	double dotp = lon.dot(X);
	double Psi = acos(dotp);
	auto dir = lon.cross(X);

	proj.row(0) = X;	// Project to orbital plane frame
	proj.row(1) = Y;
	proj.row(2) = Z;
	
	Iprj.col(0) = X;	// Project from orbital plane to original
	Iprj.col(1) = Y;
	Iprj.col(2) = Z;

	if (Nfib>1) {
				// Quadrant geometry
	  int sgn1 = 1, sgn2 = 1;
	  if (Theta > 0.5*M_PI) sgn1 = -1;
	  if (dir(2) < 0.0)     sgn2 = -1;
	  
				// Sanity check for debugging
	  if (true) {
	    double norm = 0.0;
	    if (sgn1*sgn2==1)
	      norm = (proj - return_euler(Phi, Theta,  Psi, 0)).norm();
	    else
	      norm = (proj - return_euler(Phi, Theta,  -Psi, 0)).norm();
	    
	    if (fabs(norm) > 1.0e-10) {
	      std::cout << "Theta=" << Theta
			<< " dir(2)=" << dir(2) << std::endl;
	    }
	  }

	  if (sgn1*sgn2==1) {
	    projM = return_euler(Phi, -Theta,  Psi, 0);
	    IprjM = return_euler(Phi, -Theta,  Psi, 1);
	  } else {
	    projM = return_euler(Phi, -Theta,  -Psi, 0);
	    IprjM = return_euler(Phi, -Theta,  -Psi, 1);
	  }	  
	}
      }
    }

    // Orbital plane loop
    //
    for (int q=0; q<Nrepl; q++) {

      double RR=0.0;
      for (int i=1; i<=3; i++) {
	RR += ps[i]*ps[i];
	TT += 0.5*mass*ps[0]*ps[3+i]*ps[3+i];
      }
      RR  =  sqrt(RR);
      if (RR>=rmin) {
	VC += -mass*ps[0]*rmodel->get_mass(RR)/RR;
	WW +=  0.5*mass*ps[0]*rmodel->get_pot(RR);
      }

      if (zeropos or zerovel) {
	ps[0] *= mass;
	PS.push_back(ps);
	zz[0] += ps[0];
	if (zeropos) for (int i=1; i<3; i++) zz[i] -= ps[0]*ps[i];
	if (zerovel) for (int i=4; i<7; i++) zz[i] -= ps[0]*ps[i];
      }
      else {
	out << std::setw(20) << mass * ps[0];
	for (int i=1; i<=6; i++) out << std::setw(20) << ps[i]+ps0[i];
	
	if (NI) {
	  for (int n=0; n<NI; n++) out << std::setw(10) << 0;
	}
	if (ND) {
	  for (int n=0; n<ND; n++) out << std::setw(20) << 0.0;
	}

	out << std::endl;
      }

      // Rotate particle in orbital plane
      if (Nrepl>1) {

	Eigen::Matrix3d invt;
	double Q = dq*(q/Nfib+1), phi=0.0, cost=0.0;
	int nfib = q % Nfib;

	if (Nfib>1) {
	  phi  = 2.0*M_PI * nfib / goldenRatio;
	  cost = 1.0 - 2.0 * nfib/Nfib;
	  cost = (cost < -1.0 ? -1.0 : (cost > 1.0 ? 1.0 : cost));
	  invt = return_euler(phi, acos(cost),  0, 1);
	}

	// Update rotation matrix
	rot(0, 0) = rot(1, 1) = cos(Q);
	rot(0, 1) = -sin(Q);
	rot(1, 0) =  sin(Q);

	// Full rotation matrix
	Eigen::Matrix3d trans;
	if (Nfib>1)  trans = invt * rot * proj;
	else         trans = Iprj * rot * proj;

	// Rotate position and velocity
	pos1 =  trans * pos;
	vel1 = -trans * vel;

	// Reset the phase-space vector
	//
	ps[0] = mfac;
	
	ps[1] = pos1(0);
	ps[2] = pos1(1);
	ps[3] = pos1(2);

	ps[4] = vel1(0);
	ps[5] = vel1(1);
	ps[6] = vel1(2);
      }
    }

    if (myid==0 and !((n+1)%NREPORT)) cout << '\r' << (n+1)*numprocs*Nrepl
					   << flush;
  }

  if (zeropos or zerovel) {

    if (zz[0] > 0.0) {
      for (int i=1; i<7; i++) ps0[i] += zz[i]/zz[0];
    }

    for (auto ps : PS) {
      out << std::setw(20) << ps[0];
      for (int i=1; i<=6; i++) out << std::setw(20) << ps[i]+ps0[i];

      if (NI) {
	for (int n=0; n<NI; n++) out << std::setw(10) << 0;
      }
      if (ND) {
	for (int n=0; n<ND; n++) out << std::setw(20) << 0.0;
      }

      out << std::endl;
    }
  }

  if (myid) {

    MPI_Reduce(&count, 0, 1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&negms, 0, 1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TT,    0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&WW,    0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&VC,    0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  } else {

    MPI_Reduce(MPI_IN_PLACE, &count, 1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &negms, 1, MPI_INT,    MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &TT,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &WW,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &VC,    1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    std::cout << std::endl
	      << std::setw(60) << std::setfill('-') << '-' << std::endl
	      << std::setfill(' ')
	      << std::setw(20) << "States rejected: " << count      << std::endl
	      << std::setw(20) << "Negative masses: " << negms      << std::endl
	      << std::setw(60) << setfill('-') << '-' << std::endl
	      << setfill(' ')
	      << std::setw(20) << "KE="               << TT         << std::endl
	      << std::setw(20) << "PE="               << WW         << std::endl
	      << std::setw(20) << "VC="               << VC         << std::endl
	      << std::setw(20) << "Ratio (-2T/W)="    << -2.0*TT/WW << std::endl
	      << std::setw(20) << "Ratio (-2T/C)="    << -2.0*TT/VC << std::endl
	      << std::setw(60) << std::setfill('-') << '-' << std::endl
	      << std::setfill(' ');
  }

  out.close();
  MPI_Barrier(MPI_COMM_WORLD);

  // Make the final phase-space file and clean up
  //
  if (myid==0) {
    std::ostringstream sout;
    sout << "cat " << OUTPS << ".* > " << OUTPS;
    sout << "; rm " << OUTPS << ".*";
    int ret = system(sout.str().c_str());
  }

  MPI_Finalize();

  return 0;
}
