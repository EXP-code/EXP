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

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <memory>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <fftw3.h>

#include <numerical.H>
#include <massmodel.H>
#include <model3d.H>
#include <GenPoly.H>
#include <isothermal.H>
#include <hernquist.H>
#include <EllipForce.H>
#include <localmpi.H>
#include <fpetrap.h>

namespace po = boost::program_options;

// Global variables

int
main(int argc, char **argv)
{
  int HMODEL, N, NUMDF, NUMR, NUMJ, NUME, NUMG, NREPORT, SEED, ITMAX, NUMMODEL, RNUM;
  int  DIVERGE, DIVERGE2, LOGSCALE, LINEAR, NUMINT, NI, ND;
  double DIVERGE_RFAC, DIVERGE_RFAC2, NN, MM, RA, RMODMIN, RMOD, EPS;
  double X0, Y0, Z0, U0, V0, W0, TOLE;
  double Emin0, Emax0, Kmin0, Kmax0, RBAR, MBAR, BRATIO, CRATIO, SMOOTH;
  bool LOGR, PSP, MULTI, ELIMIT, VERBOSE, GRIDPOT, MODELS, EBAR;
  std::string INFILE, MMFILE, OUTFILE, OUTPS, config;
  
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
  const char *DESC = "Generate single-mass or multi-mass spherical ICs.\n\nOptions";
  
  po::options_description desc(DESC);
  desc.add_options()
    ("help,h",
     "Print this help message")
    ("conf,c",        po::value<string>(&config),
     "Write template options file with current and all default values")
    ("input,f",       po::value<string>(&config),
     "Parameter configuration file")
    ("HMODEL",        po::value<int>(&HMODEL)->default_value(0),
     "Halo type (0=file)")
    ("N",             po::value<int>(&N)->default_value(1000000),
     "Number of bodies")
    ("NUMDF",         po::value<int>(&NUMDF)->default_value(800),
     "Number of points in energy grid for Eddington inversion")
    ("NUMR",          po::value<int>(&NUMR)->default_value(400),
     "Number of points in radial grid for halo model")
    ("NUMJ",          po::value<int>(&NUMJ)->default_value(400),
     "Number of points in kappa grid for Eddington inversion")
    ("NUME",          po::value<int>(&NUME)->default_value(400),
     "Number of energy points in PS generation table")
    ("NUMG",          po::value<int>(&NUMG)->default_value(800),
     "Number of pointsin mass table for model realization")
    ("NREPORT",       po::value<int>(&NREPORT)->default_value(1000),
     "Report after generating NREPORT points")
    ("SEED",          po::value<int>(&SEED)->default_value(11),
     "Initial seed for random number generator")
    ("ITMAX",         po::value<int>(&ITMAX)->default_value(10000),
     "Maximum number of interations for acceptance-rejection method")
    ("NUMMODEL",      po::value<int>(&NUMMODEL)->default_value(500),
     "Number of points for GeneralizedPolytrope")
    ("RNUM",          po::value<int>(&RNUM)->default_value(10000),
     "Number of radial points for interally computed mass model")
    ("DIVERGE",       po::value<int>(&DIVERGE)->default_value(0),
     "Use inner cusp extrapolation on real mass-density model")
    ("DIVERGE_RFAC",  po::value<double>(&DIVERGE_RFAC)->default_value(1.5),
     "Inner cusp slope")
    ("DIVERGE2",      po::value<int>(&DIVERGE2)->default_value(0),
     "Use inner cusp extrapolation on (pseudo) number-density model")
    ("DIVERGE_RFAC2", po::value<double>(&DIVERGE_RFAC2)->default_value(1.5),
     "Inner cusp slope")
    ("LOGR",          po::value<bool>(&LOGR)->default_value(false),
     "Use logarithmic scale for internal radial grid")
    ("LOGSCALE",      po::value<int>(&LOGSCALE)->default_value(0),
     "Use logarithmic scale for SphericalModelTable")
    ("LINEAR",        po::value<int>(&LINEAR)->default_value(0),
     "Use logarithmic scale for SphericalModelTable")
    ("NN",            po::value<double>(&NN)->default_value(2.5),
     "First polytropic index (energy)")
    ("MM",            po::value<double>(&MM)->default_value(0.5),
     "Second polytropic index (ang. mom.)")
    ("RA",            po::value<double>(&RA)->default_value(1.0e8),
     "Anisotropy index")
    ("RMODMIN",       po::value<double>(&RMODMIN)->default_value(1.0e-2),
     "Inner radius for Istothermal and Hernquist model")
    ("RMOD",          po::value<double>(&RMOD)->default_value(100.0),
     "Outer radius for Isothermal and Hernquist model")
    ("EPS",           po::value<double>(&EPS)->default_value(1.0e-5),
     "step size for computing polytrope")
    ("X0",            po::value<double>(&X0)->default_value(0.0),
     "Phase space offset")
    ("Y0",            po::value<double>(&Y0)->default_value(0.0),
     "Phase space offset")
    ("Z0",            po::value<double>(&Z0)->default_value(0.0),
     "Phase space offset")
    ("U0",            po::value<double>(&U0)->default_value(0.0),
     "Phase space offset")
    ("V0",            po::value<double>(&V0)->default_value(0.0),
     "Phase space offset")
    ("W0",            po::value<double>(&W0)->default_value(0.0),
     "Phase space offset")
    ("TOLE",          po::value<double>(&TOLE)->default_value(1.0e-4),
     "Point generation fractional energy offset for Eddington grid")
    ("Emin0",         po::value<double>(&Emin0)->default_value(-3.0),
     "Minimum energy (if ELIMIT=true)")
    ("Emax0",         po::value<double>(&Emax0)->default_value(-1.0),
     "Maximum energy (if ELIMIT=true)")
    ("Kmin0",         po::value<double>(&Kmin0)->default_value(0.0),
     "Minimum kappa (if ELIMIT=true)")
    ("Kmax0",         po::value<double>(&Kmax0)->default_value(1.0),
     "Maximum kappa (if ELIMIT=true)")
    ("RBAR",          po::value<double>(&RBAR)->default_value(0.067),
     "Semi-major axis for bar ellipsoid")
    ("MBAR",          po::value<double>(&MBAR)->default_value(0.00103739),
     "Mass of bar ellipsoid")
    ("BRATIO",        po::value<double>(&BRATIO)->default_value(0.2),
     "axis ratio b/a")
    ("CRATIO",        po::value<double>(&CRATIO)->default_value(0.05),
     "axis ratio c/b")
    ("SMOOTH",        po::value<double>(&SMOOTH)->default_value(-1.0),
     "smooth generated profile before inversion with width SMOOTH (neg. means no smoothing")
    ("NUMINT",        po::value<int>(&NUMINT)->default_value(40),
     "Number of points for bar monopole grid")
    ("PSP",           po::value<bool>(&PSP)->default_value(true),
     "PSP ascii output")
    ("MULTI",         po::value<bool>(&MULTI)->default_value(false),
     "Use multimass")
    ("ELIMIT",        po::value<bool>(&ELIMIT)->default_value(false),
     "Limit energy and angular momentum (kappa)")
    ("VERBOSE",       po::value<bool>(&VERBOSE)->default_value(false),
     "Additional diagnostic output")
    ("GRIDPOT",       po::value<bool>(&GRIDPOT)->default_value(false),
     "Compute potential internally for multimass")
    ("MODELS",        po::value<bool>(&MODELS)->default_value(false),
     "Write out internal models")
    ("EBAR",          po::value<bool>(&EBAR)->default_value(false),
     "Add an ellipsoid bar model")
    ("INFILE",        po::value<string>(&INFILE)->default_value("infile"),
     "Mass-model file for halo")
    ("MMFILE",        po::value<string>(&MMFILE),
     "Model file for  (pseudo) number-density profile")
    ("OUTFILE",       po::value<string>(&OUTFILE)->default_value("model.out"),
     "Output file for internally generated models")
    ("OUTPS",         po::value<string>(&OUTPS)->default_value("new.bods"),
     "Output phase space file")
    ("NI",            po::value<int>(&NI)->default_value(0),
     "Number of interger attributes")
    ("ND",            po::value<int>(&ND)->default_value(0),
     "Number of double attributes")
    ;
  
  po::variables_map vm;
  
  try
    {
      po::store(po::parse_command_line(argc, argv, desc), vm);
      
      po::notify(vm); // throws on error, so do after help in case
                      // there are any problems
    }
  catch(boost::program_options::error& e)
    {
      if (myid==0) std::cout << "Option error on command line: "
			     << e.what() << std::endl;
      MPI_Finalize();
      return -1;
    }
  
  // Print help message and exit
  //
  if (vm.count("help")) {
    if (myid == 0) {
      std::cout << std::string(80, '%') << std::endl
		<< desc << std::endl << std::endl
		<< "Examples: " << std::endl
		<< "\t" << "Use parameters read from a config file in INI style"  << std::endl
		<< "\t" << argv[0] << " --input=gendisk.config"  << std::endl << std::endl
		<< "\t" << "Generate a template config file in INI style from current defaults"  << std::endl
		<< "\t" << argv[0] << " --conf=template.config" << std::endl << std::endl
		<< "\t" << "Override a single parameter in a config file from the command line"  << std::endl
		<< "\t" << argv[0] << " --conf=template.config" << std::endl << std::endl;
    }
    MPI_Finalize();
    return 0;
  }
  
  // Write template config file in INI style and exit
  //
  if (vm.count("conf")) {
    // Do not overwrite existing config file
    //
    if (boost::filesystem::exists(config)) {
      if (myid == 0)
	std::cerr << argv[0] << ": config file <" << config
		  << "> exists, will not overwrite" << std::endl;
      MPI_Finalize();
      return 0;
    }

    // Write template file
    //
    if (myid==0) {
      std::ofstream out(config);

      if (out) {
	// Iterate map and print out key--value pairs and description
	//
	for (const auto& it : vm) {
				// Don't write this parameter
	  if (it.first.find("conf")==0) continue;

	  out << std::setw(20) << std::left << it.first << " = ";
	  auto& value = it.second.value();
	  if (auto v = boost::any_cast<uint32_t>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<int>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<unsigned>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<float>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<double>(&value))
	    out << std::setw(32) << std::left << *v;
	  else if (auto v = boost::any_cast<bool>(&value))
	    out << std::setw(32) << std::left << std::boolalpha << *v;
	  else if (auto v = boost::any_cast<std::string>(&value))
	    out << std::setw(32) << std::left << *v;
	  else
	    out << "error";


	  //                               NO approximations -----+
	  // Add description as a comment                         |
	  //                                                      V
	  const po::option_description& rec = desc.find(it.first, false);
	  out << " # " << rec.description() << std::endl;
	}
      } else {
	if (myid==0)
	  std::cerr << argv[0] << ": error opening template config file <"
		    << config << ">" << std::endl;
      }
    }
    MPI_Finalize();
    return 0;
  }

  // Read parameters fron the config file
  //
  if (vm.count("input")) {
    try {
      std::ifstream in(config);
      po::store(po::parse_config_file(in, desc), vm);
      po::notify(vm);
    } catch (po::error& e) {
      if (myid==0) std::cout << "Option error in configuration file: "
			     << e.what() << std::endl;
      MPI_Finalize();
      return 0;
    }
  }
  
  // Prepare output streams
  //
  std::ostringstream sout;
  sout << OUTPS << "." << myid;
  std::ofstream out(sout.str());
  if (!out) {
    std::cerr << "Couldn't open <" << sout.str() << "> for output" << std::endl;
    exit (-1);
  }
  
  out.precision(11);
  out.setf(ios::scientific, ios::floatfield);
  
  // Begin integration
  //
  AxiSymModPtr hmodel;
  SphModTblPtr htmodel, htmodel2, mmmodel;
  SphModMultPtr multi;
  
  SphericalModelTable::even = 0;
  SphericalModelTable::logscale = 1;
  
  AxiSymModel::numr     = NUMR;
  AxiSymModel::numj     = NUMJ;
  AxiSymModel::gen_N    = NUMG;
  AxiSymModel::gen_E    = NUME;
  AxiSymModel::gen_tolE = TOLE;
  
  if (HMODEL>=0) {
    
    // Halo model
    switch (HMODEL) {
    case file:
      SphericalModelTable::even = 0;
      SphericalModelTable::logscale = LOGSCALE;
      SphericalModelTable::linear   = LINEAR;
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
    
    if (myid==0) cout << "-----------" << endl;
    auto ellip = std::make_shared<EllipForce>
      (RBAR, RBAR*BRATIO, RBAR*BRATIO*CRATIO, MBAR, NUMINT, 200);
    
    if (SMOOTH>0.0) {
      
      double a[NUMR], b[NUMR], c[NUMR];
      fftw_complex A[NUMR], B[NUMR], C[NUMR];
      fftw_plan pa, pb, pinv;
      
      pa   = fftw_plan_dft_r2c_1d(NUMR, a, A, FFTW_ESTIMATE);
      pb   = fftw_plan_dft_r2c_1d(NUMR, b, B, FFTW_ESTIMATE);
      pinv = fftw_plan_dft_c2r_1d(NUMR, C, c, FFTW_ESTIMATE);
      
      double xmin=rmin, xmax=rmax;
      
      double x, z, dx = (30*SMOOTH + xmax - xmin)/NUMR;
      double y, dk = (30.0*SMOOTH + xmax - xmin)/NUMR;
      double scale = dk / NUMR;
      int indx;
      
      ofstream tin, tout;
      
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
	  tin << setw(10) << indx
	      << setw(20) << x
	      << setw(20) << a[indx]
	      << setw(20) << z
	      << setw(20) << b[indx]
	      << endl;
      }
      
      fftw_execute(pa);
      fftw_execute(pb);
      
      for (int i=0; i<NUMR; i++) {
	C[i][0] = (A[i][0] * B[i][0] - A[i][1] * B[i][1]) * scale;
	C[i][1] = (A[i][0] * B[i][1] + A[i][1] * B[i][0]) * scale;
      }
      
      // inverse transform to get c, the convolution of a and b;
      //            this has the side effect of overwriting C
      fftw_execute(pinv);
      
      // ...
      
      double mbar=-1.0;
      vector<double> rr(NUMR), mm(NUMR);
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
	  tout << setw(20) << x
	       << setw(20) << y
	       << setw(20) << mm[i]
	       << endl;
      }
      
      fftw_destroy_plan(pa);
      fftw_destroy_plan(pb);
      fftw_destroy_plan(pinv);
      
      Linear1d mass1(rr, mm);
      
      for (int i=0; i<RNUM; i++) MS[i] = mass1.eval(r2[i]);
      
    } else {
      for (int i=0; i<RNUM; i++) MS[i] = ellip->getMass(r2[i]);
    }
    
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
	outmod << "# Internal model with ellip" << endl;
	outmod << RNUM << endl;
	for (int i=0; i<RNUM; i++) {
	  if (LOGR)
	    r = rmin*exp(dr*i);
	  else
	    r = rmin + dr*i;
	  
	  outmod << setw(20) << r2[i]
		 << setw(20) << d2[i]
		 << setw(20) << m2[i]
		 << setw(20) << p2[i]
		 << endl;
	}
      } else {
	cerr << "Error opening <" << mname << endl;
      }
    }
    
  }
  
  double mass = hmodel->get_mass(hmodel->get_max_radius())/N;
  AxiSymModPtr rmodel = hmodel;
  
  if (MULTI) {
    
    SphericalModelTable::even = 0;
    SphericalModelTable::logscale = LOGSCALE;
    SphericalModelTable::linear   = LINEAR;
    
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
	cout << "-----------" << endl
	     << "Mass (computed)=" << mt2[RNUM-1] << endl;
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
	outmod << "# Internal multimass model" << endl;
	outmod << RNUM << endl;
	for (int i=0; i<RNUM; i++) {
	  if (LOGR)
	    r = rmin*exp(dr*i);
	  else
	    r = rmin + dr*i;
	  
	  outmod << setw(20) << r
		 << setw(20) << mmmodel->get_density(r)
		 << setw(20) << mmmodel->get_mass(r)
		 << setw(20) << mmmodel->get_pot(r)
		 << endl;
	}
      } else {
	cerr << "Error opening <" << mname << endl;
      }
    }
    
    // Generate the multimass model
    //
    if (EBAR)
      multi = std::make_shared<SphericalModelMulti>(htmodel2, mmmodel);
    else
      multi = std::make_shared<SphericalModelMulti>(htmodel,  mmmodel);
    
    rmodel = multi;
  }
  
  rmodel->set_seed(SEED+myid);
  rmodel->set_itmax(ITMAX);
  
  
  if (MODELS and myid==0) {
    std::ofstream outmod(OUTFILE);
    if (outmod) {
      outmod.setf(ios::scientific);
      outmod.precision(11);
      outmod << "# Internal model" << endl;
      outmod << RNUM << endl;
      for (int i=0; i<RNUM; i++) {
	if (LOGR)
	  r = rmin*exp(dr*i);
	else
	  r = rmin + dr*i;
	
	outmod << setw(20) << r
	       << setw(20) << rmodel->get_density(r)
	       << setw(20) << rmodel->get_mass(r)
	       << setw(20) << rmodel->get_pot(r)
	       << endl;
      }
    } else {
      cerr << "Error opening <" << OUTFILE << endl;
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
    if (PSP) out << setw(8) << N
		 << setw(6) << NI << setw(6) << ND << endl;
    else     out << setw(8) << N << setw(15) << 0.0 << endl;
  }
  
  // Diagnostic variables
  //
  int count=0, negms=0;
  double TT=0.0, WW=0.0, VC=0.0, RR;
  
  if (myid==0)
    cout << "-----------" << endl
	 << "Body count:" << endl;
  
  int npernode = N/numprocs;
  int beg = myid*npernode;
  int end = beg + npernode;
  
  if (myid==numprocs - 1) end = N;
  
  for (int n=beg; n<end; n++) {
    
    do {
      if (ELIMIT)
	ps = rmodel->gen_point(Emin0, Emax0, Kmin0, Kmax0, ierr);
      else		
	ps = rmodel->gen_point(ierr);
      if (ierr) count++;
    } while (ierr);
    
    out << setw(20) << mass * ps[0];
    
    if (ps[0] <= 0.0) negms++;
    
    RR=0.0;
    for (int i=1; i<=3; i++) {
      RR += ps[i]*ps[i];
      TT += 0.5*mass*ps[0]*ps[3+i]*ps[3+i];
    }
    VC += -mass*ps[0]*rmodel->get_mass(sqrt(RR))/sqrt(RR);
    WW += 0.5*mass*ps[0]*rmodel->get_pot(sqrt(RR));
    
    for (int i=1; i<=6; i++) out << setw(20) << ps[i]+ps0[i];
    
    if (PSP) {
      if (NI) {
	for (int n=0; n<NI; n++) out << setw(10) << 0;
      }
      if (ND) {
	for (int n=0; n<ND; n++) out << setw(20) << 0.0;
      }
    }
    
    out << endl;
    
    if (myid==0 and !((n+1)%NREPORT)) cout << '\r' << (n+1)*numprocs << flush;
    
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
    
    cout << endl << "-----------" << endl << endl;
    cout << setw(20) << "States rejected: " << count      << endl;
    cout << setw(20) << "Negative masses: " << negms      << endl;
    cout << setw(60) << setfill('-') << '-' << endl       << setfill(' ');
    cout << setw(20) << "KE="               << TT         << endl;
    cout << setw(20) << "PE="               << WW         << endl;
    cout << setw(20) << "VC="               << VC         << endl;
    cout << setw(20) << "Ratio (-2T/W)="    << -2.0*TT/WW << endl;
    cout << setw(20) << "Ratio (-2T/C)="    << -2.0*TT/VC << endl;
    std::cout << setw(60) << std::setfill('-') << '-' << std::endl << std::setfill(' ');
  }
  
  MPI_Finalize();
  
  return 0;
}
