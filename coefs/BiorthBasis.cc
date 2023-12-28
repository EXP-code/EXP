#include <algorithm>

#include <YamlCheck.H>
#include <EXPException.H>
#include <BiorthBasis.H>
#include <DiskModels.H>
#include <exputils.H>
#include <gaussQ.H>

#ifdef HAVE_FE_ENABLE
#include <cfenv>
#endif

namespace BasisClasses
{
  const std::set<std::string>
  SphericalSL::valid_keys = {
    "rmapping",
    "cmap",
    "Lmax",
    "dof",
    "npca",
    "npca0",
    "pcavar",
    "pcadiag",
    "pcavtk",
    "subsamp",
    "hexp",
    "snr",
    "samplesz",
    "vtkfreq",
    "tksmooth",
    "tkcum",
    "tk_type",
    "nmax",
    "modelname",
    "scale",
    "rmin",
    "rmax",
    "numr",
    "nums",
    "N1",
    "N2",
    "NO_L0",
    "NO_L1",
    "EVEN_L",
    "EVEN_M",
    "M0_ONLY",
    "NOISE",
    "noiseN",
    "noise_model_file",
    "seedN",
    "ssfrac",
    "playback",
    "coefCompute",
    "coefMaster",
    "diverge",
    "dfac",
    "dtime",
    "logr",
    "plummer",
    "self_consistent",
    "cachename"
  };

  std::vector<std::string> BiorthBasis::getFieldLabels(const Coord ctype)
  {
    // Field labels (force field labels added below)
    //
    std::vector<std::string> labels =
      {"potl", "potl m>0", "potl m=0", 
       "dens", "dens m>0", "dens m=0"};

    if (ctype == Coord::Cylindrical) {
      labels.push_back("rad force");
      labels.push_back("ver force");
      labels.push_back("azi force");
    } else if (ctype == Coord::Cartesian) {
      labels.push_back("x force");
      labels.push_back("y force");
      labels.push_back("z force");
    } else if (ctype == Coord::None) {
      // No forces
    } else {
      labels.push_back("rad force");
      labels.push_back("mer force");
      labels.push_back("azi force");
    }

    return labels;
  }

  SphericalSL::SphericalSL(const YAML::Node& CONF) : BiorthBasis(CONF)
  {
    initialize();
  }

  SphericalSL::SphericalSL(const std::string& confstr) : BiorthBasis(confstr)
  {
    initialize();
  }

  void SphericalSL::initialize()
  {

    // Assign some defaults
    //
    cmap       = 1;
    lmax       = 6;
    nmax       = 18;
    model_file = "SLGridSph.model";
    
    // Check for unmatched keys
    //
    auto unmatched = YamlCheck(conf, valid_keys);
    if (unmatched.size())
      throw YamlConfigError("Basis::Basis::SphericalSL", "parameter", unmatched, __FILE__, __LINE__);
    
    // Default cachename, empty by default
    //
    std::string cachename;

    // Assign values from YAML
    //
    double rmap = 1.0;

    try {
      if (conf["cmap"])      cmap       = conf["cmap"].as<int>();
      if (conf["Lmax"])      lmax       = conf["Lmax"].as<int>();
      if (conf["nmax"])      nmax       = conf["nmax"].as<int>();
      if (conf["modelname"]) model_file = conf["modelname"].as<std::string>();
      
      if (conf["rmapping"]) 
	rmap = conf["rmapping"].as<double>();
      else
	rmap = 1.0;
      
      if (conf["scale"]) 
	scale = conf["scale"].as<double>();
      else
	scale = 1.0;
      
      if (conf["rmin"]) 
	rmin = conf["rmin"].as<double>();
      else
	rmin = 0.0;
      
      if (conf["rmax"]) 
	rmax = conf["rmax"].as<double>();
      else
	rmax = std::numeric_limits<double>::max();
      
      if (conf["numr"])
	numr = conf["numr"].as<int>();
      else
	numr = 800;
      
      N1 = 0;
      N2 = std::numeric_limits<int>::max();
      NO_L0 = NO_L1 = EVEN_L = EVEN_M = M0_only = false;
      
      if (conf["N1"]   )     N1        = conf["N1"].as<bool>();
      if (conf["N2"]   )     N2        = conf["N2"].as<bool>();
      if (conf["NO_L0"])     NO_L0     = conf["NO_L0"].as<bool>();
      if (conf["NO_L1"])     NO_L1     = conf["NO_L1"].as<bool>();
      if (conf["EVEN_L"])    EVEN_L    = conf["EVEN_L"].as<bool>();
      if (conf["EVEN_M"])    EVEN_M    = conf["EVEN_M"].as<bool>();
      if (conf["M0_ONLY"])   M0_only   = conf["M0_ONLY"].as<bool>();
      if (conf["cachename"]) cachename = conf["cachename"].as<std::string>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameter stanza for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("SphericalSL: error parsing YAML");
    }
    
    // Check for non-null cache file name.  This must be specified
    // to prevent recomputation and unexpected behavior.
    //
    if (cachename.size() == 0) {
      throw std::runtime_error
	("SphericalSL requires a specified cachename in your YAML config\n"
	 "for consistency with previous invocations and existing coefficient\n"
	 "sets.  Please add explicitly add 'cachename: name' to your config\n"
	 "with new 'name' for creating a basis or an existing 'name' for\n"
	 "reading a previously generated basis cache\n");
    }

    // Set MPI flag in SLGridSph from MPI_Initialized
    SLGridSph::mpi = use_mpi ? 1 : 0;
    
    // Instantiate to get min/max radius from the model
    mod = std::make_shared<SphericalModelTable>(model_file);
    
    // Set rmin to a sane value if not specified
    if (not conf["rmin"] or rmin < mod->get_min_radius()) 
      rmin = mod->get_min_radius();

    // Set rmax to a sane value if not specified
    if (not conf["rmax"] or rmax > mod->get_max_radius()) 
      rmax = mod->get_max_radius()*0.99;
    
    // Finally, make the Sturm-Lioville basis...
    sl = std::make_shared<SLGridSph>
      (model_file, lmax, nmax, numr, rmin, rmax, true, cmap, rmap,
       0, 1, cachename);
    
    // Test basis for consistency
    orthoTest(orthoCheck(std::max<int>(nmax*50, 200)), classname(), harmonic());

    // Number of possible threads
    int nthrds = omp_get_max_threads();
    
    potd.resize(nthrds);
    dpot.resize(nthrds);
    dpt2.resize(nthrds);
    dend.resize(nthrds);
    
    legs  .resize(nthrds);
    dlegs .resize(nthrds);
    d2legs.resize(nthrds);

    for (auto & v : potd) v.resize(lmax+1, nmax);
    for (auto & v : dpot) v.resize(lmax+1, nmax);
    for (auto & v : dpt2) v.resize(lmax+1, nmax);
    for (auto & v : dend) v.resize(lmax+1, nmax);

    for (auto & v : legs  ) v.resize(lmax+1, lmax+1);
    for (auto & v : dlegs ) v.resize(lmax+1, lmax+1);
    for (auto & v : d2legs) v.resize(lmax+1, lmax+1);

    expcoef.resize((lmax+1)*(lmax+1), nmax);
    expcoef.setZero();
      
    work.resize(nmax);
      
    factorial.resize(lmax+1, lmax+1);
      
    for (int l=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	factorial(l, m) = sqrt( (0.5*l+0.25)/M_PI * 
				exp(lgamma(1.0+l-m) - lgamma(1.0+l+m)) );
	if (m != 0) factorial(l, m) *= M_SQRT2;
      }
    }

    used = 0;

    // Set spherical coordindates
    //
    coordinates = Coord::Spherical;
  }
  
  void SphericalSL::reset_coefs(void)
  {
    if (expcoef.rows()>0 && expcoef.cols()>0) expcoef.setZero();
    totalMass = 0.0;
    used = 0;
  }
  
  
  void SphericalSL::load_coefs(CoefClasses::CoefStrPtr coef, double time)
  {
    CoefClasses::SphStruct* cf = dynamic_cast<CoefClasses::SphStruct*>(coef.get());

    cf->lmax   = lmax;
    cf->nmax   = nmax;
    cf->scale  = scale;
    cf->time   = time;
    cf->normed = true;

    // Angular storage dimension
    int ldim = (lmax+1)*(lmax+2)/2;

    // Allocate the coefficient storage
    cf->store.resize((lmax+1)*(lmax+2)/2, nmax);

    // Make the coefficient map
    cf->coefs = std::make_shared<CoefClasses::SphStruct::coefType>
      (cf->store.data(),ldim, nmax);

    for (int l=0, L0=0, L1=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	for (int n=0; n<nmax; n++) {
	  if (m==0)
	    (*cf->coefs)(L0, n) = {expcoef(L1, n), 0.0};
	  else
	    (*cf->coefs)(L0, n) = {expcoef(L1, n), expcoef(L1+1, n)};
	}
	L0 += 1;
	if (m==0) L1 += 1;
	else      L1 += 2;
      }
    }
  }

  void SphericalSL::set_coefs(CoefClasses::CoefStrPtr coef)
  {
    // Sanity check on derived class type
    //
    if (typeid(*coef) != typeid(CoefClasses::SphStruct))
      throw std::runtime_error("SphericalSL::set_coefs: you must pass a CoefClasses::SphStruct");

    // Sanity check on dimensionality
    //
    {
      auto & cf = *dynamic_cast<CoefClasses::SphStruct*>(coef.get());

      int rows = (*cf.coefs).rows();
      int cols = (*cf.coefs).cols();
      int rexp = (lmax+1)*(lmax+2)/2;
      if (rows != rexp or cols != nmax) {
	std::ostringstream sout;
	sout << "SphericalSL::set_coefs: the basis has (lmax, nmax)=("
	     << lmax << ", " << nmax
	     << ") and the dimensions must be (rows, cols)=("
	     << rexp << ", " << nmax
	     << "). The coef structure has (rows, cols)=("
	     << rows << ", " << cols << ")";
	  
	throw std::runtime_error(sout.str());
      }
    }
    
    CoefClasses::SphStruct* cf = dynamic_cast<CoefClasses::SphStruct*>(coef.get());

    // Assign internal coefficient table (doubles) from the complex struct
    //
    for (int l=0, L0=0, L1=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	for (int n=0; n<nmax; n++) {
	  if (m==0)
	    expcoef(L1,   n) = (*cf->coefs)(L0, n).real();
	  else {
	    expcoef(L1,   n) = (*cf->coefs)(L0, n).real();
	    expcoef(L1+1, n) = (*cf->coefs)(L0, n).imag();
	  }
	}
	L0 += 1;
	if (m==0) L1 += 1;
	else      L1 += 2;
      }
    }

    // Assign center if need be
    //
    if (cf->ctr.size())
      coefctr = cf->ctr;
    else
      coefctr = {0.0, 0.0, 0.0};
  }

  void SphericalSL::accumulate(double x, double y, double z, double mass)
  {
    double fac, fac1, fac2, fac4;
    double norm = -4.0*M_PI;
    const double dsmall = 1.0e-20;
    
    // Get thread id
    int tid = omp_get_thread_num();

    //======================
    // Compute coefficients 
    //======================
    
    double r2 = (x*x + y*y + z*z);
    double r = sqrt(r2) + dsmall;
    double costh = z/r;
    double phi = atan2(y,x);
    double rs = r/scale;
    
    if (r < rmin or r > rmax) return;
    
    used++;
    totalMass += mass;
    
    sl->get_pot(potd[tid], rs);
    
    legendre_R(lmax, costh, legs[tid]);
    
    // L loop
    for (int l=0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {
      
      Eigen::VectorXd workE;
      int esize = (l+1)*nmax;
      
      // M loop
      for (int m=0, moffset=0, moffE=0; m<=l; m++) {
	
	if (m==0) {
	  fac = factorial(l, m) * legs[tid](l, m);
	  for (int n=0; n<nmax; n++) {
	    fac4 = potd[tid](l, n)*fac;
	    expcoef(loffset+moffset, n) += fac4 * norm * mass;
	  }
	  
	  moffset++;
	}
	else {
	  fac  = factorial(l, m) * legs[tid](l, m);
	  fac1 = fac*cos(phi*m);
	  fac2 = fac*sin(phi*m);
	  for (int n=0; n<nmax; n++) {
	    fac4 = potd[tid](l, n);
	    expcoef(loffset+moffset  , n) += fac1 * fac4 * norm * mass;
	    expcoef(loffset+moffset+1, n) += fac2 * fac4 * norm * mass;
	  }
	  
	  moffset+=2;
	}
      }
    }
    
  }
  
  void SphericalSL::make_coefs()
  {
    if (use_mpi) {
      
      MPI_Allreduce(MPI_IN_PLACE, &used, 1, MPI_INT,
		    MPI_SUM, MPI_COMM_WORLD);
      
      for (int l=0; l<(lmax+1)*(lmax+1); l++) {
	work = expcoef.row(l);
	MPI_Allreduce(MPI_IN_PLACE, work.data(), nmax, MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
	expcoef.row(l) = work;
      }
    }
  }
  
  std::vector<double>
  SphericalSL::sph_eval(double r, double costh, double phi)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    double fac1, cosm, sinm;
    double sinth = -sqrt(fabs(1.0 - costh*costh));
    
    fac1 = factorial(0, 0);
    
    sl->get_dens (dend[tid], r/scale);
    sl->get_pot  (potd[tid], r/scale);
    sl->get_force(dpot[tid], r/scale);
    
    legendre_R(lmax, costh, legs[tid], dlegs[tid]);
    
    double den0, pot0, potr;

    if (NO_L0) {
      den0 = 0.0;
      pot0 = 0.0;
      potr = 0.0;
    } else {
      den0 = fac1 * expcoef.row(0).dot(dend[tid].row(0));
      pot0 = fac1 * expcoef.row(0).dot(potd[tid].row(0));
      potr = fac1 * expcoef.row(0).dot(dpot[tid].row(0));
    }

    double den1 = 0.0;
    double pot1 = 0.0;
    double pott = 0.0;
    double potp = 0.0;
    
    // L loop
    for (int l=1, loffset=1; l<=lmax; loffset+=(2*l+1), l++) {
      
      // Check for even l
      if (EVEN_L and l%2) continue;
      
      // No l=1
      if (NO_L1 and l==1) continue;
      
      // M loop
      for (int m=0, moffset=0; m<=l; m++) {
	
	if (M0_only and m) continue;
	if (EVEN_M and m%2) continue;
	
	fac1 = factorial(l, m);
	if (m==0) {
	  double sumR=0.0, sumP=0.0, sumD=0.0;
	  for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	    sumR += expcoef(loffset+moffset, n) * dend[tid](l, n);
	    sumP += expcoef(loffset+moffset, n) * potd[tid](l, n);
	    sumD += expcoef(loffset+moffset, n) * dpot[tid](l, n);
	  }
	  
	  den1 += fac1*legs[tid] (l, m) * sumR;
	  pot1 += fac1*legs[tid] (l, m) * sumP;
	  potr += fac1*legs[tid] (l, m) * sumD;
	  pott += fac1*dlegs[tid](l, m) * sumP;
	  
	  moffset++;
	}
	else {
	  cosm = cos(phi*m);
	  sinm = sin(phi*m);
	  
	  double sumR0=0.0, sumP0=0.0, sumD0=0.0;
	  double sumR1=0.0, sumP1=0.0, sumD1=0.0;
	  for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	    sumR0 += expcoef(loffset+moffset+0, n) * dend[tid](l, n);
	    sumP0 += expcoef(loffset+moffset+0, n) * potd[tid](l, n);
	    sumD0 += expcoef(loffset+moffset+0, n) * dpot[tid](l, n);
	    sumR1 += expcoef(loffset+moffset+1, n) * dend[tid](l, n);
	    sumP1 += expcoef(loffset+moffset+1, n) * potd[tid](l, n);
	    sumD1 += expcoef(loffset+moffset+1, n) * dpot[tid](l, n);
	  }
	  
	  den1 += fac1 * legs[tid] (l, m) *  ( sumR0*cosm + sumR1*sinm );
	  pot1 += fac1 * legs[tid] (l, m) *  ( sumP0*cosm + sumP1*sinm );
	  potr += fac1 * legs[tid] (l, m) *  ( sumD0*cosm + sumD1*sinm );
	  pott += fac1 * dlegs[tid](l, m) *  ( sumP0*cosm + sumP1*sinm );
	  potp += fac1 * legs[tid] (l, m) *  (-sumP0*sinm + sumP1*cosm ) * m;
	  
	  moffset +=2;
	}
      }
    }
    
    double densfac = 1.0/(scale*scale*scale) * 0.25/M_PI;
    double potlfac = 1.0/scale;
    
    den1 += den0;
    pot1 += pot0;

    return
      {den0 * densfac,
       den1 * densfac,
       pot0 * potlfac,
       pot1 * potlfac,
       potr * (-potlfac)/scale,
       pott * (-potlfac),
       potp * (-potlfac)};
    //         ^
    //         |
    // Return force not potential gradient
  }


  std::vector<double>
  SphericalSL::cyl_eval(double R, double z, double phi)
  {
    double r = sqrt(R*R + z*z) + 1.0e-18;
    double costh = z/r, sinth = R/r;
    
    auto v = sph_eval(r, costh, phi);
    
    double potR = v[4]*sinth + v[5]*costh;
    double potz = v[4]*costh - v[5]*sinth;

    return {v[0], v[1], v[2], v[3], potR, potz, v[6]};
  }


  // Evaluate in cartesian coordinates
  std::vector<double>
  SphericalSL::crt_eval
  (double x, double y, double z)
  {
    double R = sqrt(x*x + y*y);
    double phi = atan2(y, x);

    auto v = cyl_eval(R, z, phi);

    //tdens0, tdens, tpotl0, tpotl, tpotR, tpotz, tpotp
    
    double tpotx = v[4]*x/R - v[6]*y/R ;
    double tpoty = v[4]*y/R + v[6]*x/R ;
    
    return {v[0], v[1], v[2], v[3], tpotx, tpoty, v[5]};
  }
  

  SphericalSL::BasisArray SphericalSL::getBasis
  (double logxmin, double logxmax, int numgrid)
  {
    // Assing return storage
    BasisArray ret (lmax+1);
    for (auto & v : ret) {
      v.resize(nmax);
      for (auto & u : v) {
	u["potential"].resize(numgrid); // Potential
	u["density"  ].resize(numgrid); // Density
	u["rforce"   ].resize(numgrid); // Radial force
      }
    }

    // Radial grid spacing
    double dx = (logxmax - logxmin)/numgrid;

    // Basis storage
    Eigen::MatrixXd tabpot, tabden, tabfrc;

    for (int i=0; i<numgrid; i++) {
      sl->get_pot  (tabpot, pow(10.0, logxmin + dx*i));
      sl->get_dens (tabden, pow(10.0, logxmin + dx*i));
      sl->get_force(tabfrc, pow(10.0, logxmin + dx*i));
      for (int l=0; l<=lmax; l++) {
	for (int n=0; n<nmax;n++){
	  ret[l][n]["potential"](i) = tabpot(l, n);
	  ret[l][n]["density"  ](i) = tabden(l, n);
	  ret[l][n]["rforce"   ](i) = tabfrc(l, n);
	}
      }
    }
    
    return ret;
  }

#define MINEPS 1.0e-10
  
  void BiorthBasis::legendre_R(int lmax, double x, Eigen::MatrixXd& p)
  {
    double fact, somx2, pll, pl1, pl2;
    
    p(0, 0) = pll = 1.0;
    if (lmax > 0) {
      somx2 = sqrt( (1.0 - x)*(1.0 + x) );
      fact = 1.0;
      for (int m=1; m<=lmax; m++) {
	pll *= -fact*somx2;
	p(m, m) = pll;
	if (std::isnan(p(m, m))) {
	  std::ostringstream sout;
	  sout << "BiorthBasis::legendre_R NaN: p[" << m << "][" << m << "]: pll=" << pll;
	  throw std::runtime_error(sout.str());
	}
	fact += 2.0;
      }
    }
    
    for (int m=0; m<lmax; m++) {
      pl2 = p(m, m);
      p(m+1, m) = pl1 = x*(2*m+1)*pl2;
      for (int l=m+2; l<=lmax; l++) {
	p(l, m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
	if (std::isnan(p(l, m))) {
	  std::ostringstream sout;
	  sout << "BiorthBasis::legendre_R NaN: p[" << l << "][" << m << "]: pll=" << pll;
	  throw std::runtime_error(sout.str());
	}

	pl2 = pl1;
	pl1 = pll;
      }
    }
    
      if (std::isnan(x)) {
	std::ostringstream sout;
	sout << "BiorthBasis::legendre_R NaN: x";
	throw std::runtime_error(sout.str());
      }
    for (int l=0; l<=lmax; l++)
      for (int m=0; m<=l; m++)
	if (std::isnan(p(l, m))) {
	  std::ostringstream sout;
	  sout << "BiorthBasis::legendre_R NaN: p[" << l << "][" << m << "] lmax=" << lmax;
	  throw std::runtime_error(sout.str());
	}
  }
  
  void BiorthBasis::legendre_R(int lmax, double x, Eigen::MatrixXd& p,
			      Eigen::MatrixXd &dp)
  {
    double fact, somx2, pll, pl1, pl2;
    
    p(0, 0) = pll = 1.0;
    if (lmax > 0) {
      somx2 = sqrt( (1.0 - x)*(1.0 + x) );
      fact = 1.0;
      for (int m=1; m<=lmax; m++) {
	pll *= -fact*somx2;
	p(m, m) = pll;
	fact += 2.0;
      }
    }
    
    for (int m=0; m<lmax; m++) {
      pl2 = p(m, m);
      p(m+1, m) = pl1 = x*(2*m+1)*pl2;
      for (int l=m+2; l<=lmax; l++) {
	p(l, m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
	pl2 = pl1;
	pl1 = pll;
      }
    }
    
    if (1.0-fabs(x) < MINEPS) {
      if (x>0) x =   1.0 - MINEPS;
      else     x = -(1.0 - MINEPS);
    }
    
    somx2 = 1.0/(x*x - 1.0);
    dp(0, 0) = 0.0;
    for (int l=1; l<=lmax; l++) {
      for (int m=0; m<l; m++)
	dp(l, m) = somx2*(x*l*p(l, m) - (l+m)*p(l-1, m));
      dp(l, l) = somx2*x*l*p(l, l);
    }
  }
  
  void BiorthBasis::legendre_R(int lmax, double x, Eigen::MatrixXd &p,
			      Eigen::MatrixXd &dp, Eigen::MatrixXd& d2p)
  {
    double fact, somx2, pll, pl1, pl2;
    
    p(0, 0) = pll = 1.0;
    if (lmax > 0) {
      somx2 = sqrt( (1.0 - x)*(1.0 + x) );
      fact = 1.0;
      for (int m=1; m<=lmax; m++) {
	pll *= -fact*somx2;
	p(m, m) = pll;
	fact += 2.0;
      }
    }
    
    for (int m=0; m<lmax; m++) {
      pl2 = p(m, m);
      p(m+1, m) = pl1 = x*(2*m+1)*pl2;
      for (int l=m+2; l<=lmax; l++) {
	p(l, m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
	pl2 = pl1;
	pl1 = pll;
      }
    }
    
    if (1.0-fabs(x) < MINEPS) {
      if (x>0) x =   1.0 - MINEPS;
      else     x = -(1.0 - MINEPS);
    }
    
    somx2 = 1.0/(x*x - 1.0);
    dp(0, 0) = 0.0;
    if (lmax) {
      for (int l=1; l<=lmax; l++) {
	for (int m=0; m<l; m++)
	  dp(l, m) = somx2*(x*l*p(l, m) - (l+m)*p(l-1, m));
	dp(l, l) = somx2*x*l*p(l, l);
      }
    }
    
    for (int l=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++)
	d2p(l, m) = -somx2*(2.0*x*dp(l, m) - p(l, m)*(somx2*m*m + l*(l+1)));
    }
    
  }
  
  // Disk target types for cylindrical basis construction
  const std::map<std::string, Cylindrical::DiskType> Cylindrical::dtlookup =
    { {"constant",    DiskType::constant},
      {"gaussian",    DiskType::gaussian},
      {"mn",          DiskType::mn},
      {"exponential", DiskType::exponential},
      {"doubleexpon", DiskType::doubleexpon}
    };

  const std::set<std::string>
  Cylindrical::valid_keys = {
    "tk_type",
    "rcylmin",
    "rcylmax",
    "acyl",
    "hcyl",
    "hexp",
    "snr",
    "evcut",
    "nmaxfid",
    "lmaxfid",
    "mmax",
    "mlim",
    "nmax",
    "ncylnx",
    "ncylny",
    "ncylr",
    "ncylodd",
    "ncylrecomp",
    "npca",
    "npca0",
    "nvtk",
    "eof_file",
    "override",
    "samplesz",
    "rnum",
    "pnum",
    "tnum",
    "ashift",
    "expcond",
    "ignore",
    "deproject",
    "logr",
    "pcavar",
    "pcaeof",
    "pcavtk",
    "pcadiag",
    "subsamp",
    "try_cache",
    "density",
    "EVEN_M",
    "cmap",
    "cmapr",
    "cmapz",
    "aratio",
    "hratio",
    "dweight",
    "rwidth",
    "rfactor",
    "rtrunc",
    "rpow",
    "mtype",
    "dtype",
    "vflag",
    "self_consistent",
    "playback",
    "coefCompute",
    "coefMaster"
  };

  Cylindrical::Cylindrical(const YAML::Node& CONF) : BiorthBasis(CONF)
  {
    initialize();
  }

  Cylindrical::Cylindrical(const std::string& confstr) : BiorthBasis(confstr)
  {
    initialize();
  }


  double Cylindrical::DiskDens(double R, double z, double phi)
  {
    double ans = 0.0;

    switch (DTYPE) {

    case DiskType::constant:
      if (R < acyl && fabs(z) < hcyl)
	ans = 1.0/(2.0*hcyl*M_PI*acyl*acyl);
      break;
      
    case DiskType::gaussian:
      if (fabs(z) < hcyl)
	ans = 1.0/(2.0*hcyl*2.0*M_PI*acyl*acyl)*
	  exp(-R*R/(2.0*acyl*acyl));
      break;
      
    case DiskType::mn:
      {
	double Z2 = z*z + hcyl*hcyl;
	double Z  = sqrt(Z2);
	double Q2 = (acyl + Z)*(acyl + Z);
	ans = 0.25*hcyl*hcyl/M_PI*(acyl*R*R + (acyl + 3.0*Z)*Q2)/( pow(R*R + Q2, 2.5) * Z*Z2 );
      }
      break;
      
    case DiskType::doubleexpon:
      {
	double a1 = acyl;
	double a2 = acyl*aratio;
	double h1 = hcyl;
	double h2 = hcyl*hratio;
	double w1 = 1.0/(1.0+dweight);
	double w2 = dweight/(1.0+dweight);
	
	double f1 = cosh(z/h1);
	double f2 = cosh(z/h2);
	
	ans =
	  w1*exp(-R/a1)/(4.0*M_PI*a1*a1*h1*f1*f1) +
	  w2*exp(-R/a2)/(4.0*M_PI*a2*a2*h2*f2*f2) ;
      }
      break;
    case DiskType::exponential:
    default:
      {
	double f = cosh(z/hcyl);
	ans = exp(-R/acyl)/(4.0*M_PI*acyl*acyl*hcyl*f*f);
      }
      break;
    }
    
    if (rwidth>0.0) ans *= erf((rtrunc-R)/rwidth);
    
    return ans;
  }

  double Cylindrical::dcond(double R, double z, double phi, int M)
  {
    //
    // No shift for M==0
    //
    if (M==0) return DiskDens(R, z, phi);
    
    //
    // Fold into [-PI/M, PI/M] for M>=1
    //
    double dmult = M_PI/M, phiS;
    if (phi>M_PI)
      phiS = phi + dmult*(int)((2.0*M_PI - phi)/dmult);
    else
      phiS = phi - dmult*(int)(phi/dmult);
    
    //
    // Apply a shift along the x-axis
    //
    double x = R*cos(phiS) - ashift*acyl;
    double y = R*sin(phiS);
    
    return DiskDens(sqrt(x*x + y*y), z, atan2(y, x));
  }


  void Cylindrical::initialize()
  {
    // Assign some defaults
    //
    rcylmin     = 0.001;
    rcylmax     = 20.0;
    acyl        = 0.01;
    hcyl        = 0.002;
    hexp        = 1.0;
    nmax        = 18;
    mmax        = 6;
    mlim        = std::numeric_limits<int>::max();
    lmaxfid     = 128;
    nmaxfid     = 64;
    ncylnx      = 256;
    ncylny      = 128;
    ncylodd     = 9;
    ncylr       = 2000;
    eof_file    = ".eof_cache_file";
    Ignore      = false;
    deproject   = false;
    
    rnum        = 200;
    pnum        = 1;
    tnum        = 80;
    ashift      = 0.0;
    logarithmic = false;
    density     = true;
    EVEN_M      = false;
    cmapR       = 1;
    cmapZ       = 1;
    mtype       = "exponential";
    dtype       = "exponential";
    vflag       = 0;
    
    // Basis construction defaults
    //
    // The EmpCylSL deprojection from specified disk model (EXP or MN) -- change only if using MN
    dmodel      = "EXP"; 

    // Radial scale length ratio for disk basis construction with doubleexpon
    aratio      = 1.0; 

    // Vertical scale height ratio for disk basis construction with doubleexpon
    hratio      = 1.0;              

    // Ratio of second disk relative to the first disk for disk basis construction with double-exponential
    dweight     = 1.0;              

    // Width for erf truncation for EOF conditioning density (ignored if zero)
    rwidth      = 0.0;             

    // Fraction of scale length for shift in conditioning function
    ashift      = 0.0;              

    // Disk radial scaling factor for spherical deprojection model
    rfactor     = 1.0;  

    // Maximum disk radius for erf truncation of EOF conditioning density
    rtrunc      = 0.1;           

    // Power-law index for power-law disk profile
    ppow        = 4.0;

    // Check for unmatched keys
    //
    auto unmatched = YamlCheck(conf, valid_keys);
    if (unmatched.size())
      throw YamlConfigError("Basis::Basis::Cylindrical", "parameter", unmatched, __FILE__, __LINE__);
    
    // Assign values from YAML
    //
    try {
      // These first two should not be user settable . . . but need them for now
      //
      if (conf["rcylmin"   ])    rcylmin  = conf["rcylmin"   ].as<double>();
      if (conf["rcylmax"   ])    rcylmax  = conf["rcylmax"   ].as<double>();
      
      if (conf["acyl"      ])       acyl  = conf["acyl"      ].as<double>();
      if (conf["hcyl"      ])       hcyl  = conf["hcyl"      ].as<double>();
      if (conf["hexp"      ])       hexp  = conf["hexp"      ].as<double>();
      if (conf["lmaxfid"   ])    lmaxfid  = conf["lmaxfid"   ].as<int>();
      if (conf["nmaxfid"   ])    nmaxfid  = conf["nmaxfid"   ].as<int>();
      if (conf["nmax"      ])       nmax  = conf["nmax"      ].as<int>();
      if (conf["mmax"      ])       mmax  = conf["mmax"      ].as<int>();
      if (conf["mlim"      ])       mlim  = conf["mlim"      ].as<int>();
      if (conf["ncylnx"    ])     ncylnx  = conf["ncylnx"    ].as<int>();
      if (conf["ncylny"    ])     ncylny  = conf["ncylny"    ].as<int>();
      if (conf["ncylr"     ])      ncylr  = conf["ncylr"     ].as<int>();
      if (conf["ncylodd"   ])    ncylodd  = conf["ncylodd"   ].as<int>();
      if (conf["eof_file"  ])   eof_file  = conf["eof_file"  ].as<std::string>();
      if (conf["rnum"      ])       rnum  = conf["rnum"      ].as<int>();
      if (conf["pnum"      ])       pnum  = conf["pnum"      ].as<int>();
      if (conf["tnum"      ])       tnum  = conf["tnum"      ].as<int>();

      if (conf["ashift"    ])     ashift  = conf["ashift"    ].as<double>();
      if (conf["logr"      ]) logarithmic = conf["logr"      ].as<bool>();
      if (conf["EVEN_M"    ])     EVEN_M  = conf["EVEN_M"    ].as<bool>();
      if (conf["cmapr"     ])      cmapR  = conf["cmapr"     ].as<int>();
      if (conf["cmapz"     ])      cmapZ  = conf["cmapz"     ].as<int>();
      if (conf["ignore"    ])      Ignore = conf["ignore"    ].as<bool>();
      if (conf["deproject" ])   deproject = conf["deproject" ].as<bool>();
      if (conf["dmodel"    ])      dmodel = conf["dmodel"    ].as<bool>();

      if (conf["aratio"    ])      aratio = conf["aratio"    ].as<double>();
      if (conf["hratio"    ])      aratio = conf["hratio"    ].as<double>();
      if (conf["dweight"   ])      aratio = conf["dweight"   ].as<double>();
      if (conf["rwidth"    ])      aratio = conf["rwidth"    ].as<double>();
      if (conf["ashift"    ])      aratio = conf["ashift"    ].as<double>();
      if (conf["rfactor"   ])     rfactor = conf["rfactor"   ].as<double>();
      if (conf["rtrunc"    ])      rtrunc = conf["rtrunc"    ].as<double>();
      if (conf["pow"       ])        ppow = conf["ppow"      ].as<double>();
      if (conf["mtype"     ])       mtype = conf["mtype"     ].as<std::string>();
      if (conf["dtype"     ])       dtype = conf["dtype"     ].as<std::string>();
      if (conf["vflag"     ])       vflag = conf["vflag"     ].as<int>();

      // Deprecation warning
      if (conf["density"   ]) {
	if (myid==0)
	  std::cout << "Cylindrical: parameter 'density' is deprecated. "
		    << "The density field will be computed regardless."
		    << std::endl;
      }

    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing 'force' for Component <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << "Config node"        << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("Cylindrical: error parsing YAML");
    }
    
    // Enforce sane values for EOF integration
    //
    rnum = std::max<int>(10, rnum);
    pnum = std::max<int>(1,  pnum);
    tnum = std::max<int>(10, tnum);
    
    EmpCylSL::RMIN        = rcylmin;
    EmpCylSL::RMAX        = rcylmax;
    EmpCylSL::NUMX        = ncylnx;
    EmpCylSL::NUMY        = ncylny;
    EmpCylSL::NUMR        = ncylr;
    EmpCylSL::CMAPR       = cmapR;
    EmpCylSL::CMAPZ       = cmapZ;
    EmpCylSL::logarithmic = logarithmic;
    EmpCylSL::VFLAG       = vflag;
    
    // Check for non-null cache file name.  This must be specified
    // to prevent recomputation and unexpected behavior.
    //
    if (not conf["eof_file"]) {
      throw std::runtime_error
	("Cylindrical requires a specified 'eof_name' in your YAML config\n"
	 "for consistency with previous invocations and existing coefficient\n"
	 "sets.  Please add explicitly add 'eof_name: name' to your config\n"
	 "with new 'name' for creating a basis or an existing 'name' for\n"
	 "reading a previously generated basis cache\n");
    }

    // Make the empirical orthogonal basis instance
    //
    sl = std::make_shared<EmpCylSL>
      (nmaxfid, lmaxfid, mmax, nmax, acyl, hcyl, ncylodd, eof_file);
    
    // Set azimuthal harmonic order restriction?
    //
    if (mlim>=0)  sl->set_mlim(mlim);
    if (EVEN_M)   sl->setEven(EVEN_M);
      
    // Attempt to read EOF cache
    //
    if (sl->read_cache() == 0) {

      // Remake cylindrical basis
      //

      if (Ignore and myid==0) {
	std::cout << "---- BasisFactory: We have deprecated the 'ignore' parameter for the" << std::endl
		  << "----               Cylindrical class. If the cache file exists but does" << std::endl
		  << "----               not match the requested parameters, the old cache file" << std::endl
		  << "----               will be moved to a .bak file, a new basis will be re-" << std::endl
		  << "----               computed, and a new cache saved in its place.  Please" << std::endl
		  << "----               remove 'ignore' from your YAML configuration." << std::endl;
      }

      // Convert mtype string to lower case
      //
      std::transform(mtype.begin(), mtype.end(), mtype.begin(),
		     [](unsigned char c){ return std::tolower(c); });
      
      // Set EmpCylSL mtype.  This is the spherical function used to
      // generate the EOF basis.  If "deproject" is set, this will be
      // overriden in EmpCylSL.
      //
      EmpCylSL::mtype = EmpCylSL::Exponential;
      if (mtype.compare("exponential")==0)
	EmpCylSL::mtype = EmpCylSL::Exponential;
      else if (mtype.compare("gaussian")==0)
	EmpCylSL::mtype = EmpCylSL::Gaussian;
      else if (mtype.compare("plummer")==0)
	EmpCylSL::mtype = EmpCylSL::Plummer;
      else if (mtype.compare("power")==0) {
	EmpCylSL::mtype = EmpCylSL::Power;
	EmpCylSL::PPOW  = ppow;
      } else {
	if (myid==0) std::cout << "No EmpCylSL EmpModel named <"
			       << mtype << ">, valid types are: "
			       << "Exponential, Gaussian, Plummer, Power "
			       << "(not case sensitive)" << std::endl;
	if (use_mpi) MPI_Finalize();
	return;
      }

      // Convert dtype string to lower case
      //
      std::transform(dtype.begin(), dtype.end(), dtype.begin(),
		     [](unsigned char c){ return std::tolower(c); });

      // Set DiskType.  This is the functional form for the disk used to
      // condition the basis.
      //
      try {				// Check for map entry, will through if the 
	DTYPE = dtlookup.at(dtype);	// key is not in the map.
	
	if (myid==0)		// Report DiskType
	  std::cout << "DiskType is <" << dtype << ">" << std::endl;
      }
      catch (const std::out_of_range& err) {
	if (myid==0) {
	  std::cout << "DiskType error in configuraton file" << std::endl;
	  std::cout << "Valid options are: ";
	  for (auto v : dtlookup) std::cout << v.first << " ";
	  std::cout << std::endl;
	}
	if (use_mpi) MPI_Finalize();
	return;
      }

      // Use these user models to deproject for the EOF spherical basis
      //
      if (deproject) {
	// The scale in EmpCylSL is assumed to be 1 so we compute the
	// height relative to the length
	//
	double H = hcyl/acyl;

	// The model instance (you can add others in DiskModels.H).
	// It's MN or Exponential if not MN.
	//
	EmpCylSL::AxiDiskPtr model;
	
	if (dmodel.compare("MN")==0) // Miyamoto-Nagai
	  model = std::make_shared<MNdisk>(1.0, H);
	else			// Default to exponential
	  model = std::make_shared<Exponential>(1.0, H);

	if (rwidth>0.0) {
	  model = std::make_shared<Truncated>(rtrunc/acyl,
					      rwidth/acyl,
					      model);
	  if (myid==0)
	    std::cout << "Made truncated model with R=" << rtrunc/acyl
		      << " and W=" << rwidth/acyl << std::endl;
	}
     
	sl->create_deprojection(H, rfactor, rnum, ncylr, model);
      }
    
      // Regenerate EOF from analytic density
      //
      std::function<double(double, double, double, int)> 
	f = [this](double R, double z, double phi, int M) -> double
	{
	  return this->dcond(R, z, phi, M);
	};

      sl->generate_eof(rnum, pnum, tnum, f);
    }

    // Orthogonality sanity check
    //
    orthoTest(orthoCheck(), classname(), harmonic());

    // Set cylindrical coordindates
    //
    coordinates = Coord::Cylindrical;
  }

  
  // Evaluate in spherical coordinates
  std::vector<double> Cylindrical::sph_eval(double r, double cth, double phi)
  {
    double sth = sqrt(1.0 - cth*cth);
    double R   = r*sth;
    double z   = r*cth;
    
    double tdens0, tdens, tpotl0, tpotl, tpotR, tpotz, tpotp;
    
    sl->accumulated_eval(R, z, phi, tpotl0, tpotl, tpotR, tpotz, tpotp);
    
    tdens = sl->accumulated_dens_eval(R, z, phi, tdens0);

    double tpotr = tpotR*R/r + tpotz*z/R ;
    double tpott = tpotR*z/r - tpotz*R/r ;

    return {tdens0, tdens, tpotl0, tpotl, tpotr, tpott, tpotp};
  }
  
  // Evaluate in cartesian coordinates
  std::vector<double> Cylindrical::crt_eval(double x, double y, double z)
  {
    double R = sqrt(x*x + y*y);
    double phi = atan2(y, x);

    double tdens0, tdens, tpotl0, tpotl, tpotR, tpotz, tpotp;

    sl->accumulated_eval(R, z, phi, tpotl0, tpotl, tpotR, tpotz, tpotp);
    
    tdens = sl->accumulated_dens_eval(R, z, phi, tdens0);

    double tpotx = tpotR*x/R - tpotp*y/R ;
    double tpoty = tpotR*y/R + tpotp*x/R ;

    return {tdens0, tdens, tpotl0, tpotl, tpotx, tpoty, tpotz};
  }
  
  // Evaluate in cylindrical coordinates
  std::vector<double> Cylindrical::cyl_eval(double R, double z, double phi)
  {
    double tdens0, tdens, tpotl0, tpotl, tpotR, tpotz, tpotp;

    sl->accumulated_eval(R, z, phi, tpotl0, tpotl, tpotR, tpotz, tpotp);
    tdens = sl->accumulated_dens_eval(R, z, phi, tdens0);

    return {tdens0, tdens, tpotl0, tpotl, tpotR, tpotz, tpotp};
  }
  
  void Cylindrical::accumulate(double x, double y, double z, double mass)
  {
    double R   = sqrt(x*x + y*y);
    double phi = atan2(y, x);
    sl->accumulate(R, z, phi, mass, 0, 0);
  }
  
  void Cylindrical::reset_coefs(void)
  {
    sl->setup_accumulation();
  }
  
  void Cylindrical::load_coefs(CoefClasses::CoefStrPtr coef, double time)
  {
    CoefClasses::CylStruct* cf = dynamic_cast<CoefClasses::CylStruct*>(coef.get());

    cf->mmax   = mmax;
    cf->nmax   = nmax;
    cf->time   = time;

    Eigen::VectorXd cos1(nmax), sin1(nmax);

    cf->store((mmax+1)*nmax);
    cf->coefs = std::make_shared<CoefClasses::CylStruct::coefType>
      (cf->store.data(), mmax+1, nmax);

    for (int m=0; m<=mmax; m++) {
      sl->get_coefs(m, cos1, sin1);

      for (int n=0; n<nmax; n++) {
	(*cf->coefs)(m, n) = {cos1(n), sin1(n)};
      }
    }
  }

  void Cylindrical::set_coefs(CoefClasses::CoefStrPtr coef)
  {
    if (typeid(*coef) != typeid(CoefClasses::CylStruct))
      throw std::runtime_error("Cylindrical::set_coefs: you must pass a CoefClasses::CylStruct");

    CoefClasses::CylStruct* cf = dynamic_cast<CoefClasses::CylStruct*>(coef.get());

    for (int m=0; m<=mmax; m++) { // Set to zero on m=0 call only--------+
      sl->set_coefs(m, (*cf->coefs).row(m).real(), (*cf->coefs).row(m).imag(), m==0);
    }

    // Assign center if need be
    //
    if (cf->ctr.size())
      coefctr = cf->ctr;
    else
      coefctr = {0.0, 0.0, 0.0};
  }

  void Cylindrical::make_coefs(void)
  {
    sl->make_coefficients();
  }
  
  
  Cylindrical::BasisArray Cylindrical::getBasis
  (double xmin, double xmax, int numR, double zmin, double zmax, int numZ)
  {
    // Allocate storage
    BasisArray ret(mmax+1);
    for (auto & v : ret) {
      v.resize(nmax);
      for (auto & u : v) {
	u["potential"].resize(numR, numZ); // Potential
	u["density"  ].resize(numR, numZ); // Density
	u["rforce"   ].resize(numR, numZ); // Radial force
	u["zforce"   ].resize(numR, numZ); // Vertical force
      }
    }
    
    // Grid spacing
    double delR = (xmax - xmin)/std::max<int>(numR-1, 1);
    double delZ = (zmax - zmin)/std::max<int>(numZ-1, 1);

    // Return values
    double p, d, fr, fz, fp;

    // Now, evaluate the grid
    for (int m=0; m<=mmax; m++) {
      for (int n=0; n<nmax; n++) {
	for (int i=0; i<numR; i++) {
	  double R = xmin + delR*i;
	  for (int j=0; j<numZ; j++) {
	    double Z = zmin + delZ*j;
	    sl->get_all(m, n, R, Z, 0.0, p, d, fr, fz, fp);
	    ret[m][n]["potential"](i,j) = p;
	    ret[m][n]["density"  ](i,j) = d;
	    ret[m][n]["rforce"   ](i,j) = fr;
	    ret[m][n]["zforce"   ](i,j) = fz;
	  }
	}
      }
    }

    return ret;
  }

  const std::set<std::string>
  FlatDisk::valid_keys = {
    "nmaxfid",
    "rcylmin",
    "rcylmax",
    "acyltbl",
    "numx",
    "numy",
    "numr",
    "knots",
    "logr",
    "model",
    "biorth",
    "scale",
    "rmin",
    "rmax",
    "self_consistent",
    "NO_M0",
    "NO_M1",
    "EVEN_M",
    "M0_ONLY",
    "ssfrac",
    "playback",
    "coefMaster",
    "Lmax",
    "Mmax",
    "nmax",
    "mmax",
    "dof",
    "subsamp",
    "samplesz",
    "vtkfreq",
    "tksmooth",
    "tkcum",
    "tk_type",
    "cachename"
  };

  FlatDisk::FlatDisk(const YAML::Node& CONF) : BiorthBasis(CONF)
  {
    initialize();
  }

  FlatDisk::FlatDisk(const std::string& confstr) : BiorthBasis(confstr)
  {
    initialize();
  }

  void FlatDisk::initialize()
  {

    // Assign some defaults
    //
    cmap       = 1;
    mmax       = 6;
    nmax       = 18;
    
    // Check for unmatched keys
    //
    auto unmatched = YamlCheck(conf, valid_keys);
    if (unmatched.size())
      throw YamlConfigError("Basis::Basis::FlatDisk", "parameter", unmatched, __FILE__, __LINE__);
    
    // Default cachename, empty by default
    //
    std::string cachename;

    // Assign values from YAML
    //
    try {
      if (conf["cmap"])      cmap       = conf["cmap"].as<int>();
      if (conf["mmax"])      mmax       = conf["mmax"].as<int>();
      if (conf["nmax"])      nmax       = conf["nmax"].as<int>();
      
      if (conf["rcylmin"]) 
	rcylmin = conf["rcylmin"].as<double>();
      else
	rcylmin = 0.0;
      
      if (conf["rcylmax"]) 
	rcylmax = conf["rcylmax"].as<double>();
      else
	rcylmax = 10.0;
      
      N1 = 0;
      N2 = std::numeric_limits<int>::max();
      NO_M0 = NO_M1 = EVEN_M = M0_only = false;
      
      if (conf["N1"]   )     N1        = conf["N1"].as<bool>();
      if (conf["N2"]   )     N2        = conf["N2"].as<bool>();
      if (conf["NO_M0"])     NO_M0     = conf["NO_M0"].as<bool>();
      if (conf["NO_M1"])     NO_M1     = conf["NO_M1"].as<bool>();
      if (conf["EVEN_M"])    EVEN_M    = conf["EVEN_M"].as<bool>();
      if (conf["M0_ONLY"])   M0_only   = conf["M0_ONLY"].as<bool>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameter stanza for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("FlatDisk: error parsing YAML");
    }
    
    // Set cmapR and cmapZ defaults
    //
    if (not conf["cmapR"])   conf["cmapR"] = cmap;
    if (not conf["cmapZ"])   conf["cmapZ"] = cmap;
    
    // Set characteristic radius defaults
    //
    if (not conf["acyltbl"]) conf["acyltbl"] = 0.6;
    if (not conf["scale"])   conf["scale"]   = 1.0;


    // Check for non-null cache file name.  This must be specified
    // to prevent recomputation and unexpected behavior.
    //
    if (not conf["cachename"]) {
      throw std::runtime_error
	("FlatDisk requires a specified cachename in your YAML config\n"
	 "for consistency with previous invocations and existing coefficient\n"
	 "sets.  Please add explicitly add 'cachename: name' to your config\n"
	 "with new 'name' for creating a basis or an existing 'name' for\n"
	 "reading a previously generated basis cache\n");
    }

    // Finally, make the basis
    //
    ortho = std::make_shared<BiorthCyl>(conf);
    
    // Orthogonality sanity check
    //
    orthoTest(orthoCheck(), classname(), harmonic());

    // Get max threads
    //
    int nthrds = omp_get_max_threads();

    // Allocate memory
    //
    potd.resize(nthrds);
    potR.resize(nthrds);
    potZ.resize(nthrds);
    dend.resize(nthrds);

    for (auto & v : potd) v.resize(mmax+1, nmax);
    for (auto & v : potR) v.resize(mmax+1, nmax);
    for (auto & v : potZ) v.resize(mmax+1, nmax);
    for (auto & v : dend) v.resize(mmax+1, nmax);

    expcoef.resize(2*mmax+1, nmax);
    expcoef.setZero();
      
    work.resize(nmax);
      
    used = 0;

    // Set cylindrical coordindates
    //
    coordinates = Coord::Cylindrical;
  }
  
  void FlatDisk::reset_coefs(void)
  {
    if (expcoef.rows()>0 && expcoef.cols()>0) expcoef.setZero();
    totalMass = 0.0;
    used = 0;
  }
  
  
  void FlatDisk::load_coefs(CoefClasses::CoefStrPtr coef, double time)
  {
    CoefClasses::CylStruct* cf = dynamic_cast<CoefClasses::CylStruct*>(coef.get());

    cf->mmax   = mmax;
    cf->nmax   = nmax;
    cf->time   = time;

    cf->store((2*mmax+1)*nmax);
    cf->coefs = std::make_shared<CoefClasses::CylStruct::coefType>
      (cf->store.data(), 2*mmax+1, nmax);

    for (int m=0, m0=0; m<=mmax; m++) {
      for (int n=0; n<nmax; n++) {
	if (m==0)
	  (*cf->coefs)(m, n) = {expcoef(m0, n), 0.0};
	else
	  (*cf->coefs)(m, n) = {expcoef(m0, n), expcoef(m0+1, n)};
      }
      if (m==0) m0 += 1;
      else      m0 += 2;
    }
  }

  void FlatDisk::set_coefs(CoefClasses::CoefStrPtr coef)
  {
    // Sanity check on derived class type
    //
    if (typeid(*coef) != typeid(CoefClasses::CylStruct))
      throw std::runtime_error("FlatDisk::set_coefs: you must pass a CoefClasses::CylStruct");

    // Sanity check on dimensionality
    //
    {
      auto cc = dynamic_cast<CoefClasses::CylStruct*>(coef.get());
      auto cf = cc->coefs;
      int rows = cf->rows();
      int cols = cf->cols();
      if (rows != mmax+1 or cols != nmax) {
	std::ostringstream sout;
	sout << "FlatDisk::set_coefs: the basis has (mmax+1, nmax)=("
	     << mmax+1 << ", " << nmax
	     << "). The coef structure has (rows, cols)=("
	     << rows << ", " << cols << ")";
	  
	throw std::runtime_error(sout.str());
      }
    }
    
    CoefClasses::CylStruct* cf = dynamic_cast<CoefClasses::CylStruct*>(coef.get());
    auto & cc = *cf->coefs;

    // Assign internal coefficient table (doubles) from the complex struct
    //
    for (int m=0, m0=0; m<=mmax; m++) {
      for (int n=0; n<nmax; n++) {
	if (m==0)
	  expcoef(m0,   n) = cc(m, n).real();
	else {
	  expcoef(m0,   n) = cc(m, n).real();
	  expcoef(m0+1, n) = cc(m, n).imag();
	}
      }
      if (m==0) m0 += 1;
      else      m0 += 2;
    }

    // Assign center if need be
    //
    if (cf->ctr.size())
      coefctr = cf->ctr;
    else
      coefctr = {0.0, 0.0, 0.0};
  }

  void FlatDisk::accumulate(double x, double y, double z, double mass)
  {
    // Normalization factors
    //
    constexpr double norm0 = 2.0*M_PI * 0.5*M_2_SQRTPI/M_SQRT2;
    constexpr double norm1 = 2.0*M_PI * 0.5*M_2_SQRTPI;

    //======================
    // Compute coefficients 
    //======================
    
    double R2 = x*x + y*y;
    double R  = sqrt(R);
    
    // Get thread id
    int tid = omp_get_thread_num();

    if (R < ortho->getRtable() and fabs(z) < ortho->getRtable()) {
    
      used++;
      totalMass += mass;
    
      double phi = atan2(y, x);

      ortho->get_pot(potd[tid], R, 0.0);
    
      // M loop
      for (int m=0, moffset=0; m<=mmax; m++) {
	
	if (m==0) {
	  for (int n=0; n<nmax; n++) {
	    expcoef(moffset, n) += potd[tid](m, n)* mass * norm0;
	  }
	  
	  moffset++;
	}
	else {
	  double ccos = cos(phi*m);
	  double ssin = sin(phi*m);
	  for (int n=0; n<nmax; n++) {
	    expcoef(moffset  , n) += ccos * potd[tid](m, n) * mass * norm1;
	    expcoef(moffset+1, n) += ssin * potd[tid](m, n) * mass * norm1;
	  }
	  moffset+=2;
	}
      }
    }
    
  }
  
  void FlatDisk::make_coefs()
  {
    if (use_mpi) {
      
      MPI_Allreduce(MPI_IN_PLACE, &used, 1, MPI_INT,
		    MPI_SUM, MPI_COMM_WORLD);
      
      for (int m=0; m<2*mmax+1; m++) {
	work = expcoef.row(m);
	MPI_Allreduce(MPI_IN_PLACE, work.data(), nmax, MPI_DOUBLE,
		      MPI_SUM, MPI_COMM_WORLD);
	expcoef.row(m) = work;
      }
    }
  }
  
  std::vector<double>FlatDisk::cyl_eval(double R, double z, double phi)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Fixed values
    constexpr double norm0 = 0.5*M_2_SQRTPI/M_SQRT2;
    constexpr double norm1 = 0.5*M_2_SQRTPI;

    double den0=0, den1=0, pot0=0, pot1=0, rpot=0, zpot=0, ppot=0;

    // Off grid evaluation
    if (R>ortho->getRtable() or fabs(z)>ortho->getRtable()) {
      double r2 = R*R + z*z;
      double r  = sqrt(r2);
      pot0 = -totalMass/r;
      rpot = -totalMass*R/(r*r2 + 10.0*std::numeric_limits<double>::min());
      zpot = -totalMass*z/(r*r2 + 10.0*std::numeric_limits<double>::min());
      
      return {den0, den1, pot0, pot1, rpot, zpot, ppot};
    }

    // Get the basis fields
    ortho->get_dens   (dend[tid],  R, z);
    ortho->get_pot    (potd[tid],  R, z);
    ortho->get_rforce (potR[tid],  R, z);
    ortho->get_zforce (potZ[tid],  R, z);
    
    // m loop
    for (int m=0, moffset=0; m<=mmax; m++) {
      
      if (m==0 and NO_M0)        continue;
      if (m==1 and NO_M1)        continue;
      if (EVEN_M and m/2*2 != m) continue;
      if (m>0  and M0_only)      break;

      if (m==0) {
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  den0 += expcoef(0, n) * dend[tid](0, n) * norm0;
	  pot0 += expcoef(0, n) * potd[tid](0, n) * norm0;
	  rpot += expcoef(0, n) * potR[tid](0, n) * norm0;
	  zpot += expcoef(0, n) * potZ[tid](0, n) * norm0;
	}
	
	moffset++;
      }  else {
	double cosm = cos(phi*m), sinm = sin(phi*m);
	double vc, vs;

	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * dend[tid](m, n);
	  vs += expcoef(moffset+1, n) * dend[tid](m, n);
	}
	
	den1 += (vc*cosm + vs*sinm) * norm1;
      
	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potd[tid](m, n);
	  vs += expcoef(moffset+1, n) * potd[tid](m, n);
	}
	
	pot1 += ( vc*cosm + vs*sinm) * norm1;
	ppot += (-vc*sinm + vs*cosm) * m * norm1;

	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potR[tid](m, n);
	  vs += expcoef(moffset+1, n) * potR[tid](m, n);
	}

	rpot += (vc*cosm + vs*sinm) * norm1;
	
	vc = vs = 0.0;
	for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	  vc += expcoef(moffset+0, n) * potZ[tid](m, n);
	  vs += expcoef(moffset+1, n) * potZ[tid](m, n);
	}

	zpot += (vc*cosm + vs*sinm) * norm1;

	moffset +=2;
      }
    }

    den1 += den0;
    pot1 += pot0;

    den0 *= -1.0;
    den1 *= -1.0;
    pot0 *= -1.0;
    pot1 *= -1.0;
    rpot *= -1.0;
    zpot *= -1.0;
    ppot *= -1.0;

    return {den0, den1, pot0, pot1, rpot, zpot, ppot};
  }


  std::vector<double> FlatDisk::sph_eval(double r, double costh, double phi)
  {
    // Cylindrical coords
    //
    double sinth = sqrt(fabs(1.0 - costh*costh));
    double R = r*sinth, z = r*costh, potR, potz;

    auto v = cyl_eval(R, z, phi);
  
    // Spherical force element converstion
    //
    double potr = potR*sinth + potz*costh;
    double pott = potR*costh - potz*sinth;

    return {v[0], v[1], v[2], v[3], potr, pott, v[6]};
  }

  std::vector<double> FlatDisk::crt_eval(double x, double y, double z)
  {
    // Cylindrical coords from Cartesian
    //
    double R = sqrt(x*x + y*y) + 1.0e-18;
    double phi = atan2(y, x);

    auto v = cyl_eval(R, z, phi);

    double potx = v[4]*x/R - v[6]*y/R;
    double poty = v[4]*y/R + v[6]*x/R;

    return {v[0], v[1], v[2], v[3], potx, poty, v[5]};
  }

  std::vector<Eigen::MatrixXd> FlatDisk::orthoCheck()
  {
    return ortho->orthoCheck();
  }
  

  FlatDisk::BasisArray FlatDisk::getBasis
  (double logxmin, double logxmax, int numgrid)
  {
    // Assing return storage
    BasisArray ret(mmax+1);
    for (auto & v : ret) {
      v.resize(nmax);
      for (auto & u : v) {
	u["potential"].resize(numgrid); // Potential
	u["density"  ].resize(numgrid); // Density
	u["rforce"   ].resize(numgrid); // Radial force
      }
    }

    // Radial grid spacing
    double dx = (logxmax - logxmin)/numgrid;

    // Basis storage
    Eigen::MatrixXd tabpot, tabden, tabrfc;

    // Evaluate on the plane
    for (int i=0; i<numgrid; i++) {
      ortho->get_pot   (tabpot, pow(10.0, logxmin + dx*i), 0.0);
      ortho->get_dens  (tabden, pow(10.0, logxmin + dx*i), 0.0);
      ortho->get_rforce(tabrfc, pow(10.0, logxmin + dx*i), 0.0);
      for (int m=0; m<=mmax; m++) {
	for (int n=0; n<nmax; n++){
	  ret[m][n]["potential"](i) = tabpot(m, n);
	  ret[m][n]["density"  ](i) = tabden(m, n);
	  ret[m][n]["rforce"   ](i) = tabrfc(m, n);
	}
      }
    }
    
    return ret;
  }

  const std::set<std::string>
  Cube::valid_keys = {
    "nminx",
    "nminy",
    "nminz",
    "nmaxx",
    "nmaxy",
    "nmaxz",
    "knots",
    "verbose",
    "check",
    "method"
  };

  Cube::Cube(const YAML::Node& CONF) : BiorthBasis(CONF)
  {
    initialize();
  }

  Cube::Cube(const std::string& confstr) : BiorthBasis(confstr)
  {
    initialize();
  }

  void Cube::initialize()
  {
    nminx = std::numeric_limits<int>::max();
    nminy = std::numeric_limits<int>::max();
    nminz = std::numeric_limits<int>::max();

    nmaxx = 6;
    nmaxy = 6;
    nmaxz = 6;

    knots = 40;

    // Check orthogonality (false by default because of its long
    // runtime and very low utility)
    //
    bool check = false;

    // Check for unmatched keys
    //
    auto unmatched = YamlCheck(conf, valid_keys);
    if (unmatched.size())
      throw YamlConfigError("Basis::Basis::Cube", "parameter", unmatched, __FILE__, __LINE__);
    
    // Default cachename, empty by default
    //
    std::string cachename;

    // Assign values from YAML
    //
    try {
      if (conf["nminx"])      nminx = conf["nminx"].as<int>();
      if (conf["nminy"])      nminy = conf["nminy"].as<int>();
      if (conf["nminz"])      nminz = conf["nminz"].as<int>();
      
      if (conf["nmaxx"])      nmaxx = conf["nmaxx"].as<int>();
      if (conf["nmaxy"])      nmaxy = conf["nmaxy"].as<int>();
      if (conf["nmaxz"])      nmaxz = conf["nmaxz"].as<int>();
      
      if (conf["knots"])      knots = conf["knots"].as<int>();

      if (conf["check"])      check = conf["check"].as<bool>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing parameter stanza for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("Cube: error parsing YAML");
    }
    
    // Finally, make the basis
    //
    ortho = std::make_shared<BiorthCube>(conf);
    
    // Orthogonality sanity check
    //
    if (check) {
      auto orth = orthoCheck();
      std::vector<Eigen::MatrixXd> mod(1);
      mod[0].resize(orth.rows(), orth.cols());
      for (int i=0; i<mod[0].size(); i++)
	mod[0].data()[i] = sqrt(std::abs(orth.data()[i]));
      
      orthoTest(mod, classname(), harmonic());
    }

    // Get max threads
    //
    int nthrds = omp_get_max_threads();

    expcoef.resize(2*nmaxx+1, 2*nmaxy+1, 2*nmaxz+1);
    expcoef.setZero();
      
    used = 0;

    // Set cartesian coordindates
    //
    coordinates = Coord::Cartesian;
  }
  
  void Cube::reset_coefs(void)
  {
    expcoef.setZero();
    totalMass = 0.0;
    used = 0;
  }
  
  
  void Cube::load_coefs(CoefClasses::CoefStrPtr coef, double time)
  {
    auto cf = dynamic_cast<CoefClasses::CubeStruct*>(coef.get());

    cf->nmaxx   = nmaxx;
    cf->nmaxy   = nmaxy;
    cf->nmaxz   = nmaxz;
    cf->time    = time;

    cf->allocate();

    *cf->coefs = expcoef;
  }

  void Cube::set_coefs(CoefClasses::CoefStrPtr coef)
  {
    // Sanity check on derived class type
    //
    if (typeid(*coef) != typeid(CoefClasses::CubeStruct))
      throw std::runtime_error("Cube::set_coefs: you must pass a CoefClasses::CubeStruct");

    // Sanity check on dimensionality
    //
    {
      auto cc = dynamic_cast<CoefClasses::CubeStruct*>(coef.get());
      auto d  = cc->coefs->dimensions();
      if (d[0] != 2*nmaxx+1 or d[1] != 2*nmaxy+1 or d[2] != 2*nmaxz+1) {
	std::ostringstream sout;
	sout << "Cube::set_coefs: the basis has (2*nmaxx+1, 2*nmaxy+1, 2*nmaxz+1)=("
	     << 2*nmaxx+1 << ", " 
	     << 2*nmaxy+1 << ", " 
	     << 2*nmaxz+1
	     << "). The coef structure has dimension=("
	     << d[0] << ", " << d[1] << ", " << d[2] << ")";
	  
	throw std::runtime_error(sout.str());
      }
    }
    
    auto cf = dynamic_cast<CoefClasses::CubeStruct*>(coef.get());
    expcoef = *cf->coefs;

    coefctr = {0.0, 0.0, 0.0};
  }

  void Cube::accumulate(double x, double y, double z, double mass)
  {
    // Truncate to cube with sides in [0,1]
    if (x<0.0)
      x += std::floor(-x) + 1.0;
    else
      x -= std::floor( x);
    
    if (y<0.0)
      y += std::floor(-y) + 1.0;
    else
      y -= std::floor( y);
    
    if (z<0.0)
      z += std::floor(-z) + 1.0;
    else
      z -= std::floor( z);
    
    
    // Recursion multipliers
    Eigen::Vector3cd step
      {std::exp(-kfac*x), std::exp(-kfac*y), std::exp(-kfac*z)};
    
    // Initial values for recursion
    Eigen::Vector3cd init
      {std::exp(-kfac*(x*nmaxx)),
       std::exp(-kfac*(y*nmaxy)),
       std::exp(-kfac*(z*nmaxz))};
    
    Eigen::Vector3cd curr(init);
    for (int ix=0; ix<=2*nmaxx; ix++, curr(0)*=step(0)) {
      curr(1) = init(1);
      for (int iy=0; iy<=2*nmaxy; iy++, curr(1)*=step(1)) {
	curr(2) = init(2);
	for (int iz=0; iz<=2*nmaxz; iz++, curr(2)*=step(2)) {
	  
	  // Compute wavenumber; recall that the coefficients are
	  // stored as: -nmax,-nmax+1,...,0,...,nmax-1,nmax
	  //
	  int ii = ix-nmaxx;
	  int jj = iy-nmaxy;
	  int kk = iz-nmaxz;

	  // Normalization
	  double norm = 1.0/sqrt(M_PI*(ii*ii + jj*jj + kk*kk));;

	  expcoef(ix, iy, iz) += - mass * curr(0)*curr(1)*curr(2) * norm;
	}
      }
    }
  }
  
  void Cube::make_coefs()
  {
    if (use_mpi) {
      
      MPI_Allreduce(MPI_IN_PLACE, &used, 1, MPI_INT,
		    MPI_SUM, MPI_COMM_WORLD);
      
      MPI_Allreduce(MPI_IN_PLACE, expcoef.data(), expcoef.size(), MPI_DOUBLE_COMPLEX,
		    MPI_SUM, MPI_COMM_WORLD);
    }
  }
  
  std::vector<double> Cube::crt_eval(double x, double y, double z)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Position vector
    Eigen::Vector3d pos {x, y, z};

    // Get the basis fields
    double den1 = ortho->get_dens(expcoef, pos).real();
    double pot1 = ortho->get_pot (expcoef, pos).real();

    auto frc = ortho->get_force(expcoef, pos);
    
    double frcx = -frc(0).real();
    double frcy = -frc(1).real();
    double frcz = -frc(2).real();

    return {0, den1, 0, pot1, frcx, frcy, frcz};
  }

  std::vector<double> Cube::cyl_eval(double R, double z, double phi)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Cartesian from Cylindrical coordinates
    double x = R*cos(phi), y = R*sin(phi);

    // Position vector
    Eigen::Vector3d pos {x, y, z};

    // Get the basis fields
    double den1 = ortho->get_dens(expcoef, pos).real();
    double pot1 = ortho->get_pot (expcoef, pos).real();

    auto frc = ortho->get_force(expcoef, pos);
    
    double frcx = frc(0).real(), frcy = frc(1).real(), frcz = frc(2).real();

    double potR =  frcx*cos(phi) + frcy*sin(phi);
    double potp = -frcx*sin(phi) + frcy*cos(phi);
    double potz =  frcz;

    potR *= -1;
    potp *= -1;
    potz *= -1;

    return {0, den1, 0, pot1, potR, potz, potp};
  }

  std::vector<double> Cube::sph_eval(double r, double costh, double phi)
  {
    // Get thread id
    int tid = omp_get_thread_num();

    // Spherical from Cylindrical coordinates
    double sinth = sqrt(fabs(1.0 - costh*costh));
    double x = r*cos(phi)*sinth, y = r*sin(phi)*sinth, z = r*costh;

    // Position vector
    Eigen::Vector3d pos {x, y, z};

    // Get the basis fields
    double den1 = ortho->get_dens(expcoef, pos).real();
    double pot1 = ortho->get_pot (expcoef, pos).real();

    auto frc = ortho->get_force(expcoef, pos);
    
    double frcx = frc(0).real();
    double frcy = frc(1).real();
    double frcz = frc(2).real();

    double potr =  frcx*cos(phi)*sinth + frcy*sin(phi)*sinth + frcz*costh;
    double pott =  frcx*cos(phi)*costh + frcy*sin(phi)*costh - frcz*sinth;
    double potp = -frcx*sin(phi)       + frcy*cos(phi);

    potr *= -1;
    pott *= -1;
    potp *= -1;
    
    return {0, den1, 0, pot1, potr, pott, potp};
  }

  Eigen::MatrixXcd Cube::orthoCheck()
  {
    return ortho->orthoCheck();
  }
  
  // Generate coeffients from a particle reader
  CoefClasses::CoefStrPtr BiorthBasis::createFromReader
  (PR::PRptr reader, std::vector<double> ctr)
  {
    CoefClasses::CoefStrPtr coef;

    if (name.compare("sphereSL") == 0)
      coef = std::make_shared<CoefClasses::SphStruct>();
    else if (name.compare("cylinder") == 0)
      coef = std::make_shared<CoefClasses::CylStruct>();
    else if (name.compare("flatdisk") == 0)
      coef = std::make_shared<CoefClasses::CylStruct>();
    else {
      std::ostringstream sout;
      sout << "Basis::createCoefficients: basis <" << name << "> not recognized"
	   << std::endl;
      throw std::runtime_error(sout.str());
    }
      
    // Is center non-zero?
    //
    bool addCenter = false;
    for (auto v : ctr) {
      if (v != 0.0) addCenter = true;
    }

    // Add the expansion center metadata
    //
    if (addCenter) coef->ctr = ctr;

    std::vector<double> pp(3), vv(3);

    reset_coefs();
    for (auto p=reader->firstParticle(); p!=0; p=reader->nextParticle()) {

      bool use = false;
      
      if (ftor) {
	pp.assign(p->pos, p->pos+3);
	vv.assign(p->vel, p->vel+3);
	use = ftor(p->mass, pp, vv, p->indx);
      } else {
	use = true;
      }

      if (use) accumulate(p->pos[0]-ctr[0],
			  p->pos[1]-ctr[1],
			  p->pos[2]-ctr[2],
			  p->mass);
    }
    make_coefs();
    load_coefs(coef, reader->CurrentTime());
    return coef;
  }

  // Generate coefficients from a phase-space table
  void BiorthBasis::initFromArray(std::vector<double> ctr)
  {
    if (name.compare("sphereSL") == 0)
      coefret = std::make_shared<CoefClasses::SphStruct>();
    else if (name.compare("cylinder") == 0)
      coefret = std::make_shared<CoefClasses::CylStruct>();
    else if (name.compare("flatdisk") == 0)
      coefret = std::make_shared<CoefClasses::CylStruct>();
    else {
      std::ostringstream sout;
      sout << "Basis::createCoefficients: basis <" << name << "> not recognized"
	   << std::endl;
      throw std::runtime_error(sout.str());
    }
      
    // Is center non-zero?
    //
    bool addCenter = false;
    for (auto v : ctr) {
      if (v != 0.0) addCenter = true;
    }

    // Add the expansion center metadata
    //
    if (addCenter) coefret->ctr = ctr;
    
    // Register the center
    //
    coefctr = ctr;

    // Clean up for accumulation
    //
    reset_coefs();
    coefindx = 0;
  }

  // Accumulate coefficient contributions from arrays
  void BiorthBasis::addFromArray(Eigen::VectorXd& m, RowMatrixXd& p, bool roundrobin)
  {
    // Sanity check: is coefficient instance created?  This is not
    // foolproof.  It is really up the user to make sure that a call
    // to initFromArray() comes first.
    //
    if (not coefret) {
      std::string msg =
	"Basis::addFromArray: you must initialize coefficient accumulation "
	"with a call to Basis::initFromArray()";
      throw std::runtime_error(msg);
    }

    std::vector<double> p1(3), v1(3, 0);

    if (p.rows() < 10 and p.cols() > p.rows()) {
      std::cout << "Basis::addFromArray: interpreting your "
		<< p.rows() << "X" << p.cols() << " input array as "
		<< p.cols() << "X" << p.rows() << "." << std::endl;

      if (p.rows()<3) {
	std::ostringstream msg;
	msg << "Basis::addFromArray: you must pass a position array with at "
	  "least three rows for x, y, z.  Yours has " << p.rows() << ".";
	throw std::runtime_error(msg.str());
      }

      for (int n=0; n<p.cols(); n++) {

	if (n % numprocs==myid or not roundrobin) {

	  bool use = true;
	  if (ftor) {
	    for (int k=0; k<3; k++) p1[k] = p(k, n);
	    use = ftor(m(n), p1, v1, coefindx);
	  } else {
	    use = true;
	  }
	  coefindx++;
	  
	  if (use) accumulate(p(0, n)-coefctr[0],
			      p(1, n)-coefctr[1],
			      p(2, n)-coefctr[2], m(n));
	}
      }
      
    } else {

      if (p.cols()<3) {
	std::ostringstream msg;
	msg << "Basis::addFromArray: you must pass a position array with at "
	  "least three columns for x, y, z.  Yours has " << p.cols() << ".";
	throw std::runtime_error(msg.str());
      }

      for (int n=0; n<p.rows(); n++) {

	if (n % numprocs==myid or not roundrobin) {

	  bool use = true;
	  if (ftor) {
	    for (int k=0; k<3; k++) p1[k] = p(n, k);
	    use = ftor(m(n), p1, v1, coefindx);
	  } else {
	    use = true;
	  }
	  coefindx++;
	  
	  if (use) accumulate(p(n, 0)-coefctr[0],
			      p(n, 1)-coefctr[1],
			      p(n, 2)-coefctr[2], m(n));
	}
      }
    }
  }

  // Generate coefficients from the accumulated array values
  CoefClasses::CoefStrPtr BiorthBasis::makeFromArray(double time)
  {
    make_coefs();
    load_coefs(coefret, time);
    return coefret;
  }

  // Generate coefficients from a phase-space table
  //
  CoefClasses::CoefStrPtr BiorthBasis::createFromArray
  (Eigen::VectorXd& m, RowMatrixXd& p, double time, std::vector<double> ctr,
   bool roundrobin)
  {
    initFromArray(ctr);
    addFromArray(m, p, roundrobin);
    return makeFromArray(time);
  }

  // This evaluation step is performed by all derived classes
  Eigen::MatrixXd& AccelFunc::evalaccel
  (Eigen::MatrixXd& ps, Eigen::MatrixXd& accel, BasisCoef mod)
  {
    // Get Model
    //
    auto basis = std::get<0>(mod);

    // Get fields
    //
    int rows = accel.rows();
    for (int n=0; n<rows; n++) {
      auto v = basis->getFields(ps(n, 0), ps(n, 1), ps(n, 2));
      for (int k=0; k<3; k++) accel(n, k) += v[4+k];
    }

    return accel;
  }

  // This is an example of a evalcoefs() derived class.  It is the
  // responsibility of the derived-class implementer to provide a sane
  // set of coefficients using Basis::set_coefs for each
  // component. Although not needed here, the best way of identifying
  // the component might be to use the getName() member of coefs,
  // e.g. 'std::string name = std::get<mod>(1)->getName();'
  void
  AllTimeAccel::evalcoefs(double t, BasisCoef mod)
  {
    auto basis = std::get<0>(mod);
    auto coefs = std::get<1>(mod);

    // Interpolate coefficients
    //
    auto times = coefs->Times();

    if (t<times.front() or t>times.back()) {
      std::ostringstream sout;
      sout << "Basis::OneAccel: time t=" << t << " is out of bounds: ["
	   << times.front() << ", " << times.back() << "]";
      throw std::runtime_error(sout.str());
    }
    
    auto it1 = std::lower_bound(times.begin(), times.end(), t);
    auto it2 = it1 + 1;

    if (it2 == times.end()) {
      it2--;
      it1 = it2 - 1;
    }

    double a = (*it2 - t)/(*it2 - *it1);
    double b = (t - *it1)/(*it2 - *it1);

    auto coefsA = coefs->getCoefStruct(*it1);
    auto coefsB = coefs->getCoefStruct(*it2);

    // Duplicate a coefficient instance
    //
    auto newcoef = coefsA->deepcopy();

    // Now interpolate the matrix
    //
    newcoef->time = t;

    auto & cN = newcoef->store;
    auto & cA = coefsA->store;
    auto & cB = coefsB->store;

    for (int i=0; i<newcoef->store.size(); i++)
      cN(i) = a * cA(i) + b * cB(i);

    // Interpolate center
    //
    if (coefsA->ctr.size() and coefsB->ctr.size()) {
      newcoef->ctr.resize(3);
      for (int k=0; k<3; k++)
	newcoef->ctr[k] = a * coefsA->ctr[k] + b * coefsB->ctr[k];
    }

    // Install coefficients
    //
    basis->set_coefs(newcoef);
  }

  SingleTimeAccel::SingleTimeAccel(double t, std::vector<BasisCoef> mod)
  {
    for (auto model : mod) {

      auto basis = std::get<0>(model);
      auto coefs = std::get<1>(model);

      // Interpolate coefficients
      //
      auto times = coefs->Times();

      if (t<times.front() or t>times.back()) {
	std::ostringstream sout;
	sout << "Basis::OneAccel: time t=" << t << " is out of bounds: ["
	     << times.front() << ", " << times.back() << "]";
	throw std::runtime_error(sout.str());
      }
      
      auto it1 = std::lower_bound(times.begin(), times.end(), t);
      auto it2 = it1 + 1;
      
      if (it2 == times.end()) {
	it2--;
	it1 = it2 - 1;
      }
      
      double a = (*it2 - t)/(*it2 - *it1);
      double b = (t - *it1)/(*it2 - *it1);
      
      auto coefsA = coefs->getCoefStruct(*it1);
      auto coefsB = coefs->getCoefStruct(*it2);
      
      // Duplicate a coefficient instance
      //
      auto newcoef = coefsA->deepcopy();
      
      // Now interpolate the matrix
      //
      newcoef->time = t;

      auto & cN = newcoef->store;
      auto & cA = coefsA ->store;
      auto & cB = coefsB ->store;

      for (int i=0; i<cN.size(); i++)
	cN(i) = a * cA(i) + b * cB(i);

      // Interpolate center
      //
      if (coefsA->ctr.size() and coefsB->ctr.size()) {
	newcoef->ctr.resize(3);
	for (int k=0; k<3; k++)
	  newcoef->ctr[k] = a * coefsA->ctr[k] + b * coefsB->ctr[k];
      }

      // Install coefficients
      //
      basis->set_coefs(newcoef);
    }
    // END: component model loop
  }
  
  //! Take one leap frog step; this can/should be generalized to a
  //! one-step class in the long run
  std::tuple<double, Eigen::MatrixXd>
  OneStep(double t, double h,
	  Eigen::MatrixXd& ps, Eigen::MatrixXd& accel,
	  std::vector<BasisCoef> bfe, AccelFunctor F)
  {
    int rows = ps.rows();

    // Drift 1/2
    for (int n=0; n<rows; n++) {
      for (int k=0; k<3; k++) ps(n, k) += ps(n, 3+k)*0.5*h;
    }

    // Kick
    accel.setZero();
    for (auto mod : bfe) {
      accel = F(t, ps, accel, mod);
    }
    for (int n=0; n<rows; n++) {
      for (int k=0; k<3; k++) ps(n, 3+k) += accel(n, k)*h;
    }
    
    // Drift 1/2
    for (int n=0; n<rows; n++) {
      for (int k=0; k<3; k++) ps(n, k) += ps(n, 3+k)*0.5*h;
    }


    return std::tuple<double, Eigen::MatrixXd>(t+h, ps);
  }


  std::tuple<Eigen::VectorXd, Eigen::Tensor<float, 3>>
  IntegrateOrbits
  (double tinit, double tfinal, double h,
   Eigen::MatrixXd ps, std::vector<BasisCoef> bfe, AccelFunctor F,
   int nout)
  {
    int rows = ps.rows();
    int cols = ps.cols();

    // ps should be a (n, 6) table of phase-space initial conditions
    //
    if (cols != 6) {
      std::ostringstream sout;
      sout << "IntegrateOrbits: phase space array should be n x 6 where n is "
	   << "the number of particles.  You specified " << cols << " columns";
      throw std::runtime_error(sout.str());
    }

    // Allocate the acceleration array
    //
    Eigen::MatrixXd accel(rows, 3);

    // Sanity check
    //
    if ( (tfinal - tinit)/h >
	 static_cast<double>(std::numeric_limits<int>::max()) )
      {
	std::cout << "BasicFactor::IntegrateOrbits: step size is too small or "
		  << "time interval is too large.\n";
	// Return empty data
	//
	return {Eigen::VectorXd(), Eigen::Tensor<float, 3>()};
      }
    
    // Number of steps
    //
    int numT = floor( (tfinal - tinit)/h );

    // Compute output step
    //
    nout = std::min<int>(numT, nout);
    double H = (tfinal - tinit)/nout;

    // Return data
    //
    Eigen::Tensor<float, 3> ret;

    try {
      ret.resize(rows, 6, nout);
    }
    catch (const std::bad_alloc& e) {
      std::cout << "BasicFactor::IntegrateOrbits: memory allocation failed: "
		<< e.what() << std::endl
		<< "Your requested number of orbits and time steps requires "
		<< floor(8.0*rows*6*nout/1e9)+1 << " GB free memory"
		<< std::endl;

      // Return empty data
      //
      return {Eigen::VectorXd(), Eigen::Tensor<float, 3>()};
    }

    // Time array
    //
    Eigen::VectorXd times(nout);
    
    // Do the work
    //
    times(0) = tinit;
    for (int n=0; n<rows; n++)
      for (int k=0; k<6; k++) ret(n, k, 0) = ps(n, k);

    double tnow = tinit;
    for (int s=1, cnt=1; s<numT; s++) {
      std::tie(tnow, ps) = OneStep(tnow, h, ps, accel, bfe, F);
      if (tnow >= H*cnt-h*1.0e-8) {
	times(cnt) = tnow;
	for (int n=0; n<rows; n++)
	  for (int k=0; k<6; k++) ret(n, k, cnt) = ps(n, k);
	cnt += 1;
      }
    }

    times(nout-1) = tnow;
    for (int n=0; n<rows; n++)
      for (int k=0; k<6; k++) ret(n, k, nout-1) = ps(n, k);
    
    return {times, ret};
  }

}
// END namespace BasisClasses
