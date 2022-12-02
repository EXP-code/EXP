#include <YamlCheck.H>
#include <EXPException.H>
#include <BasisFactory.H>
#include <gaussQ.H>

namespace Basis
{
  
  Basis::Basis(const YAML::Node& CONF)
  {
    // Copy the YAML config
    //
    node = CONF;

    // Complete the initialization
    //
    initialize();
  }

  Basis::Basis(const std::string& confstr)
  {
    try {
      // Read the YAML from a string
      //
      node = YAML::Load(confstr);
    }
    catch (const std::runtime_error& error) {
      std::cout << "Basis constructor: found a problem in the YAML config"
		<< std::endl;
      throw;
    }

    // Complete the initialization
    //
    initialize();
  }

  void Basis::initialize()
  {
    // Check whether MPI is initialized
    //
    int flag;
    MPI_Initialized(&flag);
    if (flag) use_mpi = true;
    else      use_mpi = false;
    
    // Fall back sanity (works for me but this needs to be fixed
    // generally)
    //
    if (use_mpi) {
      MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
      MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    }

    // Parameters for force
    //
    try {
      // First get the id name; we know it exists because it has been
      // parsed by the factory
      name = node["id"].as<std::string>();
      // Then . . . 
      conf = node["parameters"];
    }
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing force 'parameters' for <"
			     << name << ">: "
			     << error.what() << std::endl
			     << std::string(60, '-') << std::endl
			     << node                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("Basis: error parsing YAML");
    }
    
  }
  
  const std::set<std::string>
  SphericalSL::valid_keys = {
    "cmap",
    "Lmax",
    "nmax",
    "modelname",
    "scale",
    "rmin",
    "rmax",
    "numr",
    "N1",
    "N2",
    "NO_L0",
    "NO_L1",
    "EVEN_L",
    "EVEN_M",
    "M0_ONLY"
  };

  SphericalSL::SphericalSL(const YAML::Node& CONF) : Basis(CONF)
  {
    initialize();
  }

  SphericalSL::SphericalSL(const std::string& confstr) : Basis(confstr)
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
    
    // Assign values from YAML
    //
    try {
      if (conf["cmap"])      cmap       = conf["cmap"].as<int>();
      if (conf["Lmax"])      lmax       = conf["Lmax"].as<int>();
      if (conf["nmax"])      nmax       = conf["nmax"].as<int>();
      if (conf["modelname"]) model_file = conf["modelname"].as<std::string>();
      
      if (conf["scale"]) 
	rscl = conf["scale"].as<double>();
      else
	rscl = 1.0;
      
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
      
      if (conf["N1"]   )   N1      = conf["N1"].as<bool>();
      if (conf["N2"]   )   N2      = conf["N2"].as<bool>();
      if (conf["NO_L0"])   NO_L0   = conf["NO_L0"].as<bool>();
      if (conf["NO_L1"])   NO_L1   = conf["NO_L1"].as<bool>();
      if (conf["EVEN_L"])  EVEN_L  = conf["EVEN_L"].as<bool>();
      if (conf["EVEN_M"])  EVEN_M  = conf["EVEN_M"].as<bool>();
      if (conf["M0_ONLY"]) M0_only = conf["M0_ONLY"].as<bool>();
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
    
    SLGridSph::mpi = use_mpi ? 1 : 0;
    
    mod = std::make_shared<SphericalModelTable>(model_file);
    
    rmin = max<double>(mod->get_min_radius()*2.0, 
		       mod->get_max_radius()*1.0e-4);
    rmax = mod->get_max_radius()*0.99;
    
    sl = std::make_shared<SLGridSph>
      (mod, lmax, nmax, numr, rmin, rmax, true, cmap, rscl);
    
    rscl = 1.0;
    
    potd.resize(lmax+1, nmax);
    dpot.resize(lmax+1, nmax);
    dpt2.resize(lmax+1, nmax);
    dend.resize(lmax+1, nmax);
    
    legs  .resize(lmax+1, lmax+1);
    dlegs .resize(lmax+1, lmax+1);
    d2legs.resize(lmax+1, lmax+1);

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
  }
  
  void SphericalSL::reset_coefs(void)
  {
    if (expcoef.rows()>0 && expcoef.cols()>0) expcoef.setZero();
    used = 0;
  }
  
  
  void SphericalSL::load_coefs(Coefs::CoefStrPtr coef, double time)
  {
    Coefs::SphStruct* cf = dynamic_cast<Coefs::SphStruct*>(coef.get());

    cf->lmax   = lmax;
    cf->nmax   = nmax;
    cf->scale  = rscl;
    cf->time   = time;
    cf->normed = true;

    cf->coefs.resize((lmax+1)*(lmax+2)/2, nmax);
    for (int l=0, L0=0, L1=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	for (int n=0; n<nmax; n++) {
	  if (m==0)
	    cf->coefs(L0, n) = {expcoef(L1, n), 0.0};
	  else
	    cf->coefs(L0, n) = {expcoef(L1, n), expcoef(L1+1, n)};
	}
	L0 += 1;
	if (m==0) L1 += 1;
	else      L1 += 2;
      }
    }
  }

  void SphericalSL::set_coefs(Coefs::CoefStrPtr coef)
  {
    Coefs::SphStruct* cf = dynamic_cast<Coefs::SphStruct*>(coef.get());

    // Assign internal coefficient table (doubles) from the complex struct
    //
    for (int l=0, L0=0, L1=0; l<=lmax; l++) {
      for (int m=0; m<=l; m++) {
	for (int n=0; n<nmax; n++) {
	  if (m==0)
	    expcoef(L1,   n) = cf->coefs(L0, n).real();
	  else {
	    expcoef(L1,   n) = cf->coefs(L0, n).real();
	    expcoef(L1+1, n) = cf->coefs(L0, n).imag();
	  }
	}
	L0 += 1;
	if (m==0) L1 += 1;
	else      L1 += 2;
      }
    }
  }

  void SphericalSL::accumulate(double x, double y, double z, double mass)
  {
    double fac, fac1, fac2, fac4;
    double norm = -4.0*M_PI;
    const double dsmall = 1.0e-20;
    
    //======================
    // Compute coefficients 
    //======================
    
    double r2 = (x*x + y*y + z*z);
    double r = sqrt(r2) + dsmall;
    double costh = z/r;
    double phi = atan2(y,x);
    double rs = r/rscl;
    
    if (r < rmin or r > rmax) return;
    
    used++;
    totalMass += mass;
    
    sl->get_pot(potd, rs);
    
    legendre_R(lmax, costh, legs);
    
    // L loop
    for (int l=0, loffset=0; l<=lmax; loffset+=(2*l+1), l++) {
      
      Eigen::VectorXd workE;
      int esize = (l+1)*nmax;
      
      // M loop
      for (int m=0, moffset=0, moffE=0; m<=l; m++) {
	
	if (m==0) {
	  fac = factorial(l, m) * legs(l, m);
	  for (int n=0; n<nmax; n++) {
	    fac4 = potd(l, n)*fac;
	    expcoef(loffset+moffset, n) += fac4 * norm * mass;
	  }
	  
	  moffset++;
	}
	else {
	  fac  = factorial(l, m) * legs(l, m);
	  fac1 = fac*cos(phi*m);
	  fac2 = fac*sin(phi*m);
	  for (int n=0; n<nmax; n++) {
	    fac4 = potd(l, n);
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
  
  void SphericalSL::getFields
  (double x, double y, double z,
   double& tdens0, double& tpotl0, double& tdens, double& tpotl, 
   double& tpotx, double& tpoty, double& tpotz)
  {
    double r = sqrt(x*x + y*y + z*z) + 1.0e-18;
    double theta = acos(z/r);
    double phi   = atan2(y, x);
    double cth   = cos(theta), sth = sin(theta);
    double cph   = cos(phi),   sph = sin(phi);
    double tpotr, tpott, tpotp;
    
    all_eval(r, cth, phi,
	     tdens0, tdens, tpotl0, tpotl, tpotr, tpott, tpotp);
    
    tpotx = tpotr*sth*cph + tpott*cth*cph - tpotp*sph;
    tpoty = tpotr*sth*sph + tpott*cth*sph + tpotp*cph;
    tpotz = tpotr*cth     - tpott*sth;
  }
  
  void SphericalSL::all_eval(double r, double costh, double phi,
			     double& den0, double& den1,
			     double& pot0, double& pot1,
			     double& potr, double& pott, double& potp)
  {
    double fac1, cosm, sinm;
    double sinth = -sqrt(fabs(1.0 - costh*costh));
    
    fac1 = factorial(0, 0);
    
    sl->get_dens (dend, r/rscl);
    sl->get_pot  (potd, r/rscl);
    sl->get_force(dpot, r/rscl);
    
    legendre_R(lmax, costh, legs, dlegs);
    
    if (NO_L0) {
      den0 = 0.0;
      pot0 = 0.0;
      potr = 0.0;
    } else {
      den0 = fac1 * expcoef.row(0).dot(dend.row(0));
      pot0 = fac1 * expcoef.row(0).dot(potd.row(0));
      potr = fac1 * expcoef.row(0).dot(dpot.row(0));
    }
    den1 = 0.0;
    pot1 = 0.0;
    pott = 0.0;
    potp = 0.0;
    
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
	    sumR += expcoef(loffset+moffset, n) * dend(l, n);
	    sumP += expcoef(loffset+moffset, n) * potd(l, n);
	    sumD += expcoef(loffset+moffset, n) * dpot(l, n);
	  }
	  
	  den1 += fac1*legs (l, m) * sumR;
	  pot1 += fac1*legs (l, m) * sumP;
	  potr += fac1*legs (l, m) * sumD;
	  pott += fac1*dlegs(l, m) * sumP;
	  
	  moffset++;
	}
	else {
	  cosm = cos(phi*m);
	  sinm = sin(phi*m);
	  
	  double sumR0=0.0, sumP0=0.0, sumD0=0.0;
	  double sumR1=0.0, sumP1=0.0, sumD1=0.0;
	  for (int n=std::max<int>(0, N1); n<=std::min<int>(nmax-1, N2); n++) {
	    sumR0 += expcoef(loffset+moffset+0, n) * dend(l, n);
	    sumP0 += expcoef(loffset+moffset+0, n) * potd(l, n);
	    sumD0 += expcoef(loffset+moffset+0, n) * dpot(l, n);
	    sumR1 += expcoef(loffset+moffset+1, n) * dend(l, n);
	    sumP1 += expcoef(loffset+moffset+1, n) * potd(l, n);
	    sumD1 += expcoef(loffset+moffset+1, n) * dpot(l, n);
	  }
	  
	  den1 += fac1 * legs (l, m) *  ( sumR0*cosm + sumR1*sinm );
	  pot1 += fac1 * legs (l, m) *  ( sumP0*cosm + sumP1*sinm );
	  potr += fac1 * legs (l, m) *  ( sumD0*cosm + sumD1*sinm );
	  pott += fac1 * dlegs(l, m) *  ( sumP0*cosm + sumP1*sinm );
	  potp += fac1 * legs (l, m) *  (-sumP0*sinm + sumP1*cosm ) * m;
	  
	  moffset +=2;
	}
      }
    }
    
    double densfac = 1.0/(rscl*rscl*rscl) * 0.25/M_PI;
    double potlfac = 1.0/rscl;
    
    den0  *=  densfac;
    den1  *=  densfac;
    pot0  *=  potlfac;
    pot1  *=  potlfac;
    potr  *= -potlfac/rscl;
    pott  *= -potlfac*sinth;
    potp  *= -potlfac;
    //       ^
    //       |
    //       +--- Return force not potential gradient
  }


  std::vector<Eigen::MatrixXd> SphericalSL::orthoCheck(int num)
  {
    // Gauss-Legendre knots and weights
    LegeQuad lw(num);

    // Get the scaled coordinate limits
    double ximin = sl->r_to_xi(rmin);
    double ximax = sl->r_to_xi(rmax);

    // Initialize the return matrices
    std::vector<Eigen::MatrixXd> ret(lmax+1);
    for (auto & v : ret) v.resize(nmax, nmax);

    // Do each harmonic order
    for (int L=0; L<=lmax; L++) {

      // Unroll the loop for OpenMP parallelization
#pragma omp parallel for
      for (int nn=0; nn<nmax*nmax; nn++) {
	int n1 = nn/nmax;
	int n2 = nn - n1*nmax;

	// The inner product
	double ans=0.0;
	for (int i=0; i<num; i++) {
	  
	  double x = ximin + (ximax - ximin)*lw.knot(i);
	  double r = sl->xi_to_r(x);
	  
	  ans += r*r*sl->get_pot(x, L, n1, 0)*
	    sl->get_dens(x, L, n2, 0) /
	    sl->d_xi_to_r(x) * (ximax - ximin)*lw.weight(i);
	  
	}
	// END: inner product
	    
	// Assign the matrix element
	//
	ret[L](n1, n2) = - ans;
	//               ^
	//               |
	//               +--- Switch to normed scalar product rather
	//                    that normed gravitational energy
      }
      // END: unrolled loop
    }
    // END: harmonic order loop

    return ret;
  }

  
  std::vector<std::vector<Eigen::VectorXd>> SphericalSL::getBasis
  (double logxmin, double logxmax, int numgrid)
  {
    // Assing return storage
    std::vector<std::vector<Eigen::VectorXd>> ret(lmax+1);
    for (auto & v : ret) {
      v.resize(nmax);
      for (auto & u : v) u.resize(numgrid);
    }

    // Radial grid spacing
    double dx = (logxmax - logxmin)/numgrid;

    // Basis storage
    Eigen::MatrixXd tab;

    for (int i=0; i<numgrid; i++) {
      sl->get_pot(tab, pow(10.0, logxmin + dx*i));
      for (int l=0; l<=lmax; l++) {
	for (int n=0; n<nmax;n++){
	  ret[l][n](i) = tab(l, n);
	}
      }
    }
    
    return ret;
  }

#define MINEPS 1.0e-10
  
  void Basis::legendre_R(int lmax, double x, Eigen::MatrixXd& p)
  {
    double fact, somx2, pll, pl1, pl2;
    
    p(0, 0) = pll = 1.0;
    if (lmax > 0) {
      somx2 = sqrt( (1.0 - x)*(1.0 + x) );
      fact = 1.0;
      for (int m=1; m<=lmax; m++) {
	pll *= -fact*somx2;
	p(m, m) = pll;
	if (std::isnan(p(m, m)))
	  cerr << "legendre_R: p[" << m << "][" << m << "]: pll=" << pll << "\n";
	fact += 2.0;
      }
    }
    
    for (int m=0; m<lmax; m++) {
      pl2 = p(m, m);
      p(m+1, m) = pl1 = x*(2*m+1)*pl2;
      for (int l=m+2; l<=lmax; l++) {
	p(l, m) = pll = (x*(2*l-1)*pl1-(l+m-1)*pl2)/(l-m);
	if (std::isnan(p(l, m)))
	  cerr << "legendre_R: p[" << l << "][" << m << "]: pll=" << pll << "\n";
	
	pl2 = pl1;
	pl1 = pll;
      }
    }
    
    if (std::isnan(x))
      cerr << "legendre_R: x\n";
    for (int l=0; l<=lmax; l++)
      for (int m=0; m<=l; m++)
	if (std::isnan(p(l, m)))
	  cerr << "legendre_R: p[" << l << "][" << m << "] lmax=" 
	       << lmax << "\n";
    
  }
  
  void Basis::legendre_R(int lmax, double x, Eigen::MatrixXd& p,
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
  
  void Basis::legendre_R(int lmax, double x, Eigen::MatrixXd &p,
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
  
  const std::set<std::string>
  Cylindrical::valid_keys = {
    "rcylmin",
    "rcylmax",
    "acyl",
    "hcyl",
    "hexp",
    "lmax",
    "nmax",
    "mmax",
    "mlim",
    "ncylnx",
    "ncylny",
    "ncylr",
    "ncylorder",
    "ncylodd",
    "eof_file",
    "override",
    "rnum",
    "pnum",
    "tnum",
    "ashift",
    "logr",
    "try_cache",
    "density",
    "EVEN_M",
    "cmapr",
    "cmapz"
  };

  Cylindrical::Cylindrical(const YAML::Node& CONF) : Basis(CONF)
  {
    initialize();
  }

  Cylindrical::Cylindrical(const std::string& confstr) : Basis(confstr)
  {
    initialize();
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
    lmax        = 128;
    nmax        = 64;
    mlim        = std::numeric_limits<int>::max();
    ncylnx      = 256;
    ncylny      = 128;
    ncylorder   = 18;
    ncylodd     = 9;
    eof_file    = ".eof_cache_file";
    eof_over    = false;
    
    rnum        = 100;
    pnum        = 1;
    tnum        = 40;
    ashift      = 0.025;
    logarithmic = false;
    try_cache   = true;
    density     = true;
    EVEN_M      = false;
    cmapR       = 1;
    cmapZ       = 1;
    
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
      if (conf["lmax"      ])       lmax  = conf["lmax"      ].as<int>();
      if (conf["nmax"      ])       nmax  = conf["nmax"      ].as<int>();
      if (conf["mmax"      ])       mmax  = conf["mmax"      ].as<int>();
      if (conf["mlim"      ])       mlim  = conf["mlim"      ].as<int>();
      if (conf["ncylnx"    ])     ncylnx  = conf["ncylnx"    ].as<int>();
      if (conf["ncylny"    ])     ncylny  = conf["ncylny"    ].as<int>();
      if (conf["ncylr"     ])      ncylr  = conf["ncylr"     ].as<int>();
      if (conf["ncylorder" ])  ncylorder  = conf["ncylorder" ].as<int>();
      if (conf["ncylodd"   ])    ncylodd  = conf["ncylodd"   ].as<int>();
      if (conf["eof_file"  ])   eof_file  = conf["eof_file"  ].as<std::string>();
      if (conf["override"  ])   eof_over  = conf["override"  ].as<bool>();
      
      if (conf["rnum"      ])       rnum  = conf["rnum"      ].as<int>();
      if (conf["pnum"      ])       pnum  = conf["pnum"      ].as<int>();
      if (conf["tnum"      ])       tnum  = conf["tnum"      ].as<int>();
      if (conf["ashift"    ])     ashift  = conf["ashift"    ].as<double>();
      if (conf["logr"      ]) logarithmic = conf["logr"      ].as<bool>();
      if (conf["try_cache" ])  try_cache  = conf["try_cache" ].as<bool>();
      if (conf["density"   ])    density  = conf["density"   ].as<bool>();
      if (conf["EVEN_M"    ])     EVEN_M  = conf["EVEN_M"    ].as<bool>();
      if (conf["cmapr"     ])      cmapR  = conf["cmapr"     ].as<int>();
      if (conf["cmapz"     ])      cmapZ  = conf["cmapz"     ].as<int>();
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
    
    // Default cache file name
    //
    std::string cachename = outdir + ".eof.cache." + runtag;
    
    // EOF default file name override.  Default uses runtag suffix as
    // above.  Override file must exist if explicitly specified.
    //
    if (eof_file.size()) cachename = eof_file;
    
    // For debugging; no use by force algorithm
    //
    if (density) EmpCylSL::DENS = true;
    
    
    // Make the empirical orthogonal basis instance
    //
    sl = std::make_shared<EmpCylSL>
      (nmax, lmax, mmax, ncylorder, acyl, hcyl, ncylodd, cachename);
    
    // Set azimuthal harmonic order restriction?
    //
    if (mlim>=0)  sl->set_mlim(mlim);
    if (EVEN_M)   sl->setEven(EVEN_M);
    
    // Set up the coefficient storage and read basis
    //
    sl->read_cache();
  }
  
  void Cylindrical::getFields
  (double x, double y, double z,
   double& tdens0, double& tpotl0, double& tdens, double& tpotl, 
   double& tpotx, double& tpoty, double& tpotz)
  {
    double R   = sqrt(x*x + y*y);
    double phi = atan2(y, x);
    double cph = cos(phi), sph = sin(phi);
    
    double tpotR, tpotP;
    
    sl->accumulated_eval(R, z, phi, tpotl0, tpotl, tpotR, tpotz, tpotP);
    
    tdens = sl->accumulated_dens_eval(R, z, phi, tdens0);
    
    tpotx = tpotR*cph - tpotP*sph ;
    tpoty = tpotR*sph + tpotP*cph ;
  }
  
  // Evaluate in on spherical coordinates (should this be Cartesian)
  void Cylindrical::all_eval
  (double r, double cth, double phi,
   double& tdens0, double& tdens, double& tpotl0, double& tpotl, 
   double& tpotr, double& tpott, double& tpotp)
  {
    double sth = sqrt(1.0 - cth*cth);
    double R   = r*sth;
    double z   = r*cth;
    
    double tpotR, tpotz;
    
    sl->accumulated_eval(R, z, phi, tpotl0, tpotl, tpotR, tpotz, tpotp);
    
    tdens = sl->accumulated_dens_eval(R, z, phi, tdens0);
    tdens -= tdens0;
    tpotl -= tpotl0;

    tpotr = tpotR*R/r + tpotz*z/R ;
    tpott = tpotR*z/r - tpotz*R/r ;
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
  
  void Cylindrical::load_coefs(Coefs::CoefStrPtr coef, double time)
  {
    Coefs::CylStruct* cf = dynamic_cast<Coefs::CylStruct*>(coef.get());

    cf->mmax   = mmax;
    cf->nmax   = nmax;
    cf->time   = time;

    Eigen::VectorXd cos1(nmax), sin1(nmax);

    cf->coefs.resize(mmax+1, nmax);

    for (int m=0; m<=mmax; m++) {
      sl->get_coefs(m, cos1, sin1);

      for (int n=0; n<nmax; n++) {
	cf->coefs(m, n) = {cos1(n), sin1(n)};
      }
    }
  }

  void Cylindrical::set_coefs(Coefs::CoefStrPtr coef)
  {
    Coefs::CylStruct* cf = dynamic_cast<Coefs::CylStruct*>(coef.get());

    for (int m=0; m<=mmax; m++) { // Set to zero on m=0 call only--------+
      sl->set_coefs(m, cf->coefs.row(m).real(), cf->coefs.row(m).imag(), m==0);
    }
  }

  void Cylindrical::make_coefs(void)
  {
    sl->make_coefficients();
  }
  
  
  std::vector<std::vector<Eigen::MatrixXd>> Cylindrical::getBasis
  (double xmin, double xmax, int numR, double zmin, double zmax, int numZ)
  {
    // Allocate storage
    std::vector<std::vector<Eigen::MatrixXd>> ret(mmax+1);
    for (auto & v : ret) {
      v.resize(ncylorder);
      for (auto & u : v) u.resize(numR, numZ);
    }
    
    // Grid spacing
    double delR = (xmax - xmin)/std::max<int>(numR-1, 1);
    double delZ = (zmax - zmin)/std::max<int>(numZ-1, 1);

    // Return values
    double p, d, fr, fz, fp;

    // Now, evaluate the grid
    for (int m=0; m<=mmax; m++) {
      for (int n=0; n<ncylorder; n++) {
	for (int i=0; i<numR; i++) {
	  double R = xmin + delR*i;
	  for (int j=0; j<numZ; j++) {
	    double Z = zmin + delZ*j;
	    sl->get_all(m, n, R, Z, 0.0, p, d, fr, fz, fp);
	    ret[m][n](i,j) = p;
	  }
	}
      }
    }

    return ret;
  }

  std::shared_ptr<Basis> Basis::factory_string(const std::string& conf)
  {
    YAML::Node node;
    
    try {
      // Read the YAML from a string
      //
      node = YAML::Load(conf);
    }
    catch (const std::runtime_error& error) {
      std::cout << "Basis constructor: found a problem in the YAML config"
		<< std::endl;
      throw;
    }

    // Complete the initialization
    //
    return factory_initialize(node);
  }

  std::shared_ptr<Basis> Basis::factory(const YAML::Node& conf)
  {
    return factory_initialize(conf);
  }


  std::shared_ptr<Basis> Basis::factory_initialize(const YAML::Node& conf)
  {
    std::shared_ptr<Basis> basis;
    std::string name;
    
    // Load parameters from YAML configuration node
    try {
      name = conf["id"].as<std::string>();
    } 
    catch (YAML::Exception & error) {
      if (myid==0) std::cout << "Error parsing force id in BasisFactory"
			     << std::string(60, '-') << std::endl
			     << conf                 << std::endl
			     << std::string(60, '-') << std::endl;
      
      throw std::runtime_error("BasisFactory: error parsing YAML");
    }
    
    try {
      if ( !name.compare("sphereSL") ) {
	basis = std::make_shared<SphericalSL>(conf);
      }
      else if ( !name.compare("cylinder") ) {
	basis = std::make_shared<Cylindrical>(conf);
      }
      else {
	std::string msg("I don't know about the basis named: ");
	msg += name;
	msg += ". Known types are currently 'sphereSL' and 'cylinder'";
	throw std::runtime_error(msg);
      }
    }
    catch (std::exception& e) {
      std::cout << "Error in BasisFactory constructor: " << e.what() << std::endl;
      throw;			// Rethrow the exception?
    }
    
    return basis;
  }

  // Generate coeffients from a particle reader
  Coefs::CoefStrPtr Basis::createFromReader
  (PR::PRptr reader, std::vector<double> ctr)
  {
    Coefs::CoefStrPtr coef;

    if (name.compare("sphereSL") == 0)
      coef = std::make_shared<Coefs::SphStruct>();
    else if (name.compare("cylinder") == 0)
      coef = std::make_shared<Coefs::CylStruct>();
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
  void Basis::initFromArray(std::vector<double> ctr)
  {
    if (name.compare("sphereSL") == 0)
      coefret = std::make_shared<Coefs::SphStruct>();
    else if (name.compare("cylinder") == 0)
      coefret = std::make_shared<Coefs::CylStruct>();
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
  void Basis::addFromArray(Eigen::VectorXd& m, RowMatrixXd& p)
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

    for (int n=0; n<p.rows(); n++) {
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

  // Generate coefficients from the accumulated array values
  Coefs::CoefStrPtr Basis::makeFromArray(double time)
  {
    make_coefs();
    load_coefs(coefret, time);
    return coefret;
  }

  // Generate coefficients from a phase-space table
  //
  // I'm leaving the original version ad the default until the new
  // version is tested
#if 0
  Coefs::CoefStrPtr Basis::createFromArray
  (Eigen::VectorXd& m, RowMatrixXd& p, double time, std::vector<double> ctr)
  {
    initFromArray(ctr);
    addFromArray(m, p);
    return makeFromArray(time);
  }
#else
  Coefs::CoefStrPtr Basis::createFromArray
  (Eigen::VectorXd& m, RowMatrixXd& p, double time, std::vector<double> ctr)
  {
    Coefs::CoefStrPtr coef;

    if (name.compare("sphereSL") == 0)
      coef = std::make_shared<Coefs::SphStruct>();
    else if (name.compare("cylinder") == 0)
      coef = std::make_shared<Coefs::CylStruct>();
    else {
      std::ostringstream sout;
      sout << "Basis::createFromArray: basis <" << name << "> not recognized"
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

    reset_coefs();
    std::vector<double> p1(3), v1(3, 0);
    unsigned long indx = 0;

    for (int n=0; n<p.rows(); n++) {
      bool use = true;
      if (ftor) {
	for (int k=0; k<3; k++) p1[k] = p(n, k);
	use = ftor(m(n), p1, v1, indx);
      } else {
	use = true;
      }

      if (use) accumulate(p(n, 0)-ctr[0], p(n, 1)-ctr[1], p(n, 2)-ctr[2], m(n));
    }
    make_coefs();
    load_coefs(coef, time);
    return coef;
  }
#endif

}
// END namespace Basis
