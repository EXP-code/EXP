
#include <Disk2d.H>

// Valid key list for Disk2d including BiorthCyl keys
//
std::set<std::string> Disk2d::valid_keys = {
  "acyltbl",
  "rcylmin",
  "rcylmax",
  "numr",
  "knots",
  "logr",
  "model",
  "biorth",
  "Mmax",
  "nmaxfid",
  "nmax",
  "numx",
  "numy",
  "NQDHT",
  "cmapR",
  "cmapZ",
  "scale",
  "cachename",
  "verbose",
  "mmin",
  "mlim",
  "nmin",
  "nlim",
  "EVEN_M"
};

Disk2d::Disk2d(const YAML::Node& conf) : conf(conf)
{
				// Defaults
  knots      = 40;
  numr       = 2000;
  acyltbl    = 1.0;
  rcylmin    = 0.0;
  rcylmax    = 10.0;
  logr       = false;

				// Get initialization info
  initialize();

  if (myid==0) {
    std::string sep("----    ");
    std::cout << "---- Disk2d parameters: "
	      << std::endl << sep << "mmax="        << mmax
	      << std::endl << sep << "nmax="        << nmax
	      << std::endl << sep << "acyltbl="     << acyltbl
	      << std::endl << sep << "rcylmin="     << rcylmin
	      << std::endl << sep << "rcylmax="     << rcylmax
	      << std::endl << sep << "logr="        << std::boolalpha << logr
	      << std::endl;
  }

}


void Disk2d::initialize()
{
  // Make set from the entire db
  //
  std::set<std::string> current_keys;
  for (YAML::const_iterator it=conf.begin(); it!=conf.end(); ++it) {
    current_keys.insert(it->first.as<std::string>());
  }

  // Remove matched keys
  //
  for (auto v : valid_keys) current_keys.erase(v);
  
  // Have a fit?
  //
  if (current_keys.size()) {
    std::ostringstream sout;
    sout << "Disk2d::initialize found the following invalid keys <";
    for (auto v : current_keys) sout << ' ' << v;
    sout << ">";
    throw std::runtime_error(sout.str());
  }

  // Assign values from YAML
  //
  try {
    if (conf["acyltbl"])   acyltbl    = conf["acyltbl"].as<double>();
    if (conf["rcylmin"])   rcylmin    = conf["rcylmin"].as<double>();
    if (conf["rcylmax"])   rcylmax    = conf["rcylmax"].as<double>();
    if (conf["Mmax"])      mmax       = conf["Mmax"].as<int>();
    if (conf["nmax"])      nmax       = conf["nmax"].as<int>();
    if (conf["numr"])      numr       = conf["numr"].as<int>();
    if (conf["knots"])     knots      = conf["knots"].as<int>();
    if (conf["logr"])      logr       = conf["logr"].as<bool>();
    if (conf["model"])     model      = conf["model"].as<std::string>();
    if (conf["biorth"])    biorth     = conf["biorth"].as<std::string>();
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in Disk2d: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf                 << std::endl
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  // Create the BiorthCyl instance
  ortho = std::make_shared<BiorthCyl>(conf);

  // Acculation
  accum_setup = false;
  coefs_made = false;
}

Disk2d::~Disk2d(void)
{
  // NADA
}


void Disk2d::setup_accum()
{
  if (accum_setup) return;

  // The maximum number of possible threads for this process
  nthrds = omp_get_max_threads();

  potd.resize(nthrds);
  for (auto & v : potd) v.resize(2*mmax+1, nmax);
    
  dend.resize(2*mmax+1, nmax);
  potl.resize(2*mmax+1, nmax);
  potr.resize(2*mmax+1, nmax);
  potz.resize(2*mmax+1, nmax);

  cylmass1.resize(nthrds);
  expcoef0.resize(nthrds);
  for (auto & v : expcoef0) v.resize(2*mmax+1, nmax);
    
  expcoef.resize(2*mmax+1, nmax);
  
  accum_setup = true;
}

void Disk2d::setup_accumulation()
{
  setup_accum();

  // Zero accumulation vectors
  for (auto & t : expcoef0) t.setZero();

  coefs_made = false;
}

void Disk2d::get_dpotl(double r, double z,
			 Eigen::MatrixXd& p,
			 Eigen::MatrixXd& dpr,
			 Eigen::MatrixXd& dpz, int tid)
{
  ortho->get_pot   (p,   r, z);
  ortho->get_rforce(dpr, r, z);
  ortho->get_zforce(dpz, r, z);
}

void Disk2d::get_potl(double r, double z, Eigen::MatrixXd& p, int tid)
{
  ortho->get_pot(p, r, z);
}

void Disk2d::get_dens(double r, double z, Eigen::MatrixXd& p, int tid)
{
  ortho->get_dens(p, r, z);
}

void Disk2d::get_potl_dens(double r, double z, Eigen::MatrixXd& p,
			     Eigen::MatrixXd& d, int tid)
{
  ortho->get_pot (p, r, z);
  ortho->get_dens(d, r, z);
}


void Disk2d::accumulate(double R, double phi, double mass)
{
  constexpr double norm0 = -2.0*M_SQRT2/M_2_SQRTPI;
  constexpr double norm1 = -4.0/M_2_SQRTPI;

  double r = R + 10.0*std::numeric_limits<double>::min();
      
  if (r < getRtable()) {

    int id = omp_get_thread_num();

    get_potl(r, 0.0, potd[id], id);
    
    // m loop
    //
    for (int m=0, moffset=0; m<=mmax; m++) {

      if (m==0) {
	for (int n=0; n<nmax; n++) {
	  expcoef0[id](moffset, n) += potd[id](m, n)*mass*norm0;
	}
	
	moffset++;
      }
      else {
	double fac1 = cos(phi*m);
	double fac2 = sin(phi*m);
	
	for (int n=0; n<nmax; n++) {
	  
	  double wk = potd[id](m, n)*mass*norm1;
	  
	  expcoef0[id](moffset  , n) += wk*fac1;
	  expcoef0[id](moffset+1, n) += wk*fac2;
	}
	
	moffset+=2;
      } // m!=0
      
    } // m loop
    
    cylmass1[id] += mass;
    
  } // r < rmax

}


void Disk2d::accumulate(std::vector<Particle>& part)
{
  if (coefs_made) setup_accumulation();

  std::vector<Particle>::iterator it;
#pragma omp parallel for default(none) shared(part)
  for (it=part.begin(); it < part.end(); it++) {
    double x = it->pos[0], y = it->pos[1];
    accumulate(sqrt(x*x + y*y), atan2(y, x), it->mass);
  }
}

void Disk2d::make_coefficients()
{
  // Reduce threads
  for (int i=1; i<nthrds; i++) {
    expcoef0[0] += expcoef0[i];
    cylmass1[0] += cylmass1[i];
  }

  // MPI reduce
  MPI_Allreduce(expcoef0[0].data(), expcoef.data(), expcoef.size(), MPI_DOUBLE,
		MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&cylmass1[0], &cylmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


std::tuple<double, double, double, double, double, double, double>
Disk2d::accumulated_eval(double R, double z, double phi)
{
  constexpr double norm0 = M_2_SQRTPI/(2.0*M_SQRT2);
  constexpr double norm1 = M_2_SQRTPI/2.0;

  double d0=0.0, d1=0.0, p0=0.0, p1=0.0, fr=0.0, fz=0.0, fp=0.0;

  if (R>getRtable() or fabs(z)>getRtable()) {
    
    double r2 = R*R + z*z;
    double r  = sqrt(r2);
    p0 = -cylmass/r;
    fr = -cylmass*R/(r*r2);
    fz = -cylmass*z/(r*r2);
    
  } else {

    ortho->get_dens  (dend, R, z);
    ortho->get_pot   (potl, R, z);
    ortho->get_rforce(potr, R, z);
    ortho->get_zforce(potz, R, z);
    
    // m loop
    for (int m=0, moffset=0; m<=mmax; m++) {

      if (m==0) {
	for (int n=0; n<nmax; n++) {
	  d0 += expcoef(moffset, n)*norm0 * dend(m, n);
	  p0 += expcoef(moffset, n)*norm0 * potl(m, n);
	  fr += expcoef(moffset, n)*norm0 * potr(m, n);
	  fz += expcoef(moffset, n)*norm0 * potz(m, n);
	}
	
	moffset++;

      } else {
	
	double cosm = cos(phi*m), sinm = sin(phi*m);

	for (int n=0; n<nmax; n++) {
	  d1 += ( expcoef(moffset, n)*cosm + expcoef(moffset+1, n)*sinm )*norm1*
	    dend(m, n);
	  p1 += ( expcoef(moffset, n)*cosm + expcoef(moffset+1, n)*sinm )*norm1*
	    potl(m, n);
	  fr += ( expcoef(moffset, n)*cosm + expcoef(moffset+1, n)*sinm )*norm1*
	    potr(m, n);
	  fz += ( expcoef(moffset, n)*cosm + expcoef(moffset+1, n)*sinm )*norm1*
	    potz(m, n);
	  fp += (-expcoef(moffset, n)*sinm + expcoef(moffset+1, n)*cosm )*norm1*
	    potl(m, n);
	}

	moffset += 2;
      }
    }
  }
    
  return {d0, d1, p0, p1, fr, fz, fp};
}


void Disk2d::write_coefficients(const std::string& outfile)
{
  std::ofstream out(outfile);

  if (out) {

    out << "#" << std::string(60, '-') << std::endl
	<< "# cylmass=" << cylmass << std::endl
	<< "#" << std::string(60, '-') << std::endl;

    for (int m=0, moffset=0; m<=mmax; m++) {

      if (m==0) {
	  
	for (int n=0; n<nmax; n++) {
	  out  << std::setw(4)  << m
	       << std::setw(4)  << n
	       << std::setw(18) << expcoef(moffset, n)
	       << std::endl;
	}

	moffset++;
	
      } else {
	  
	for (int n=0; n<nmax; n++) {
	  out << std::setw(4)  << m
	      << std::setw(4)  << n
	      << std::setw(18) << expcoef(moffset  , n)
	      << std::setw(18) << expcoef(moffset+1, n)
	      << std::endl;
	}
	
	moffset+=2;
      }
    }
    // END: m loop
    out << std::string(60, '-') << std::endl;
  } else {
    std::cout << "Disk2d::write_coefficients: could not open file <"
	      << outfile << ">" << std::endl;
  }
}


// Compute the m=0 coefficients for the provided density profile
void Disk2d::inner_product(std::function<double(int, double)> dens,
			   int numr, int nlim)
{
  setup_accumulation();

  LegeQuad lege(numr);

  double xmin = ortho->r_to_xi(rcylmin);
  double xmax = ortho->r_to_xi(rcylmax);
  double dx   = xmax - xmin;

  expcoef.setZero();
  
  // Normalization:
  //
  // 1. Factor of -2*pi for density in biorthonormality
  //
  // 2. A factor of 2*pi (for m=0 which integrand of 1) and pi (for
  //    m>0 which is integrand of cos(m*phi)^2 or sin(m*phi)^2) for
  //    integral over azimuthal angle
  //
  // 3. A factor of 1/sqrt(2*pi) (for m=0) and 1/sqrt(pi) (for m>0)
  //    for complex coefficient norm
  //
  constexpr double normD = -2.0*M_PI;

  for (int l=0; l<numr; l++) {

    double x   = xmin + dx*lege.knot(l);
    double r   = ortho->xi_to_r(x);
    double jac = dx*lege.weight(l)/ortho->d_xi_to_r(x);
    
    for (int m=0, moffset=0; m<=mmax; m++) {

      double norm = normD * 2.0*M_PI/sqrt(2.0*M_PI);
      if (m) norm = normD *     M_PI/sqrt(    M_PI);

      for (int n=0; n<std::min<int>(nmax, nlim); n++)
	expcoef(moffset, n) +=
	  jac * r * dens(m, r) * ortho->get_pot(r, 0.0, m, n) * norm;

      if (m) moffset += 2;
      else   moffset += 1;
    }
  }

  coefs_made = true;
}
