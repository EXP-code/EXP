#include <iostream>
#include <iomanip>
#include <memory>

#include <yaml-cpp/yaml.h>	// YAML support

#include <Progress.H>		// Progress bar

#include <DiskEval.H>

#ifdef HAVE_OMP_H
#include <omp.h>
#endif

const std::string DiskEval::cachefile = ".DiskEval.cache";

DiskEval::DiskEval
(EmpCylSL::AxiDiskPtr model, double rmin, double rmax, double ascl,
 int lmax, int numr, int nint, bool use_progress, int mmax, int nump, bool cache) :
  model(model), ascl(ascl), lmax(lmax), numr(numr)
{


  cout << "MMAX=" << mmax << "  NUMP=" << nump << "  ascl=" << ascl << endl;

  // choose a maximum m order. pass this in later.
  //int mmax = 16;

  
  if (ascl>0.0) xscl = true;
  else          xscl = false;

  // Assign grid parameters
  //
  if (xscl) {			// x = r/a(1 + r/a) scaling
    xmin = r_to_x(rmin);
    xmax = r_to_x(rmax);
    logr = false;
  } else {
    if (rmin > 0.0) {		// Log scaling
      xmin = log(rmin);
      xmax = log(rmax);
      logr = true;
    } else {			// Linear scaling
      xmin = rmin;
      xmax = rmax;
      logr = false;
    }
  }

  dx = (xmax - xmin)/(numr-1);	// Radial grid spacing

  if (cache) {
    if (read_cache()) return;
  }

  // First: compute the disk density components
  //
  rho.resize(lmax+1);
  //for (auto & v : rho) v.resize(numr, 0.0);
    for (int m=0;m<=lmax;m++) {
    rho[m].resize(lmax+1);
      for (auto & v : rho[m]) v.resize(numr, 0.0);
  }
    
  // Gauss-Legendre
  //
  std::shared_ptr<LegeQuad> lr = std::make_shared<LegeQuad>(nint);

  int numthrd = 1, tid = 0;
#ifdef HAVE_OMP_H
#pragma omp parallel
  {
    numthrd = omp_get_num_threads();
    int myid = omp_get_thread_num();
    if (myid==0)
      std::cout << "Number of threads=" << numthrd << std::endl;
  }
#endif

  std::shared_ptr<progress::progress_display> progress;
  if (use_progress) {
    std::cout << std::endl << "Begin: exact force evaluation"
	      << std::endl << "-----------------------------"
	      << std::endl;
      
    std::cout << std::endl << "Multipole density coefficients"
	      << std::endl;
    progress = std::make_shared<progress::progress_display>(numr/numthrd);
  }
  
  // Compute \rho_{lm} on radial grid
  //
#pragma omp parallel for
  for (int i=0; i<numr; i++) {
#ifdef HAVE_OMP_H
    tid = omp_get_thread_num();
#endif

    double x = xmin + dx*i, r;

    if (xscl) {
      r = x_to_r(x);
    } else {
      if (logr) r = exp(x);
      else      r = x;
    }

    // where does the phi loop go in here?
    
    // hardwire in the number of phi wedges for integration --
    //
    double dphi = 2.*M_PI/nump;
    
    for (int p=0; p<nump; p++) {
      double phi = p*dphi;
    
      for (int n=0; n<nint; n++) {
        double cosx = lr->knot(n); // Assume rho(R, phi, z) = rho(R,phi, -z);
        double R = sqrt(1.0 - cosx*cosx) * r;
        double z = cosx * r;

        double dens = (*model)(R, z, phi);
        // override the density function here with function defined above
        //double dens = return_density(R, phi, z);

        double fac = 2.0*lr->weight(n) * dphi;

	// rho[l][m][r] += Ylm(l, m, cosx) * dens(R,z) * 2.0*M_PI * 2.0*lr->weight(n)
	
        // Even terms only--------+
        //                        |
        //                        v
        for (int l=0; l<=lmax; l+=2) {
	  int ll = l/2;

	  // set up a limit on mmax
	  int maxm = std::min<int>(mmax, l);
	  
	  // Only integrate positive m
	  //             | 
	  //             | 
	  //             v
	  for (int m=0; m<=maxm; m++) {
	
	    //if ((p==0) && (n<10)) std::cout << setw(14) << p << setw(14) << n << setw(14) << dens << endl;

	    rho[ll][m][i] += Ylm(l, m, cosx) * dens * fac * cos(m*phi);

	  } // end m loop for Ylm
        } // end l loop for Ylm
      } // end radial quadrature loop
    } // end phi quadrature loop

    // Progress bar
    //
    if (progress and tid==0) ++(*progress);
  }

  // Compute Term1 and Term2 by quadrature
  //
  T1.resize(lmax+1);
    for (int m=0;m<=lmax;m++) {
    T1[m].resize(lmax+1);
      for (auto & v : T1[m]) v.resize(numr, 0.0);
  }
    
  T2.resize(lmax+1);
   for (int m=0;m<=lmax;m++) {
    T2[m].resize(lmax+1);
      for (auto & v : T2[m]) v.resize(numr, 0.0);
  }

  int numrads = 0;
  for (int l=0; l<=lmax; l+=2) numrads++;
  numrads *= numr;
  numrads /= numthrd;

  if (use_progress) {
    std::cout << std::endl << "Quadrature loop multipole expansion"
	      << std::endl;
    progress = std::make_shared<progress::progress_display>(numrads);
  }
  
  // l loop
  //
  for (int l=0; l<=lmax; l+=2) {
    int ll = l/2;

    // m loop
    int maxm = std::min<int>(mmax, l);
    for (int m=0; m<=maxm; m+=1) {
      

    // Outer r loop
    //
#pragma omp parallel for
    for (int i=0; i<numr; i++) {
#ifdef HAVE_OMP_H
      tid = omp_get_thread_num();
#endif

      double xx = xmin + dx*i, rr;

      if (xscl) {
	rr = x_to_r(xx);
      } else {
	if (logr) rr = exp(xx);
	else      rr = xx;
      }

      // Integral over 0 to r
      //
      double sum1 = 0.0, sum2 = 0.0;

      if (i>0) {
	for (int j=0; j<=i; j++) {

	  double xl = xmin + dx*(j-1);
	  double xp = xmin + dx*j;
	  double rl = xl;
	  double rp = xp;

	  if (xscl) {
	    rl = x_to_r(xl);
	    rp = x_to_r(xp);
	  } else {
	    if (logr) {
	      rl = exp(xl);
	      rp = exp(xp);
	    }
	  }
	  
	  // Trapezoidal rule for Term 1
	  //
	  if (xscl) {
	    sum1 += 0.5*(rho[ll][m][j-1]*pow(rl/rr, l+2)*dr_to_dx(xl) +
			 rho[ll][m][j  ]*pow(rp/rr, l+2)*dr_to_dx(xp)) * dx;
	  } else {
	    if (logr) {
	      sum1 += 0.5*(rho[ll][m][j-1]*pow(rl/rr, l+3) +
			   rho[ll][m][j  ]*pow(rp/rr, l+3)) * rr * dx;
	    } else {
	      sum1 += 0.5*(rho[ll][m][j-1]*pow(rl/rr, l+2) +
			   rho[ll][m][j  ]*pow(rp/rr, l+2)) * dx;
	    }
	  }
	}
      }
	
      // Integral over r to inf
      //
      for (int j=i; j<numr-1; j++) {

	double xl = xmin + dx*j;
	double xp = xmin + dx*(j+1);
	double rl = xl;
	double rp = xp;

	if (xscl) {
	  rl = x_to_r(xl);
	  rp = x_to_r(xp);
	} else {
	  if (logr) {
	    rl = exp(xl);
	    rp = exp(xp);
	  }
	}

	// Trapezoidal rule
	//
	if (xscl) {
	  sum2 += 0.5*(rho[ll][m][j  ]*pow(rl/rr, 1-l) * dr_to_dx(xl) +
		       rho[ll][m][j+1]*pow(rp/rr, 1-l) * dr_to_dx(xp)) * dx;
	} else {
	  if (logr) {
	    sum2 += 0.5*(rho[ll][m][j  ]*pow(rl/rr, 2-l) +
			 rho[ll][m][j+1]*pow(rp/rr, 2-l)) * rr * dx;
	  } else {
	    sum2 += 0.5*(rho[ll][m][j  ]*pow(rl/rr, 1-l) +
			 rho[ll][m][j+1]*pow(rp/rr, 1-l)) * dx;
	  }
	}
      }

      // Save integral values
      //
      T1[ll][m][i] = sum1;
      T2[ll][m][i] = sum2;

      // Progress bar
      //
      if (progress and tid==0) ++(*progress);
    }
    // END: outer loop over r

    }
    // END: loop over m
    
  }
  // END: loop over l
  
  if (progress and tid==0) std::cout << std::endl;


  // Test output
  //
  std::ofstream test("DiskEval.rho");
  if (test) {

    for (int i=0; i<numr; i++) {

      double r, x = xmin + dx*i;

      if (xscl) {
	r = x_to_r(x);
      } else {
	if (logr) r = exp(x);
	else      r = x;
      }

      test << std::setw(18) << r;
      for (int l=0; l<=lmax; l+=2) {
	int ll = l/2;
	for (int m=0; m<=l; m++) 
	test << std::setw(18) << rho[ll][m][i];
      }
      test << std::endl;
    }
  }

  test.close();
  test.open("DiskEval.t1");
  if (test) {
    for (int i=0; i<numr; i++) {

      double r, x = xmin + dx*i;

      if (xscl) {
	r = x_to_r(x);
      } else {
	if (logr) r = exp(x);
	else      r = x;
      }

      test << std::setw(18) << r;
      for (int l=0; l<=lmax; l+=2) {
	int ll = l/2;
	for (int m=0;m<=l;m++) 
	test << std::setw(18) << T1[ll][m][i];
      }
      test << std::endl;
    }
  }

  test.close();
  test.open("DiskEval.t2");
  if (test) {
    for (int i=0; i<numr; i++) {
      double r, x = xmin + dx*i;
      if (xscl) {
	r = x_to_r(x);
      } else {
	if (logr) r = exp(x);
	else      r = x;
      }

      test << std::setw(18) << r;
      for (int l=0; l<=lmax; l+=2) {
	int ll = l/2;
	for (int m=0; m<=l; m++)
	test << std::setw(18) << T2[ll][m][i];
      }
      test << std::endl;
    }
  }

  write_cache();
}


std::tuple<double, double, double, double> DiskEval::operator()(double R,
double z, double phi)
{
  // Get spherical coordinates
  //
  double r     = sqrt(R*R + z*z);
  double cosx  = z/r;

  // Grid interpolation values
  //
  double x;
  int i1, i2;

  if (xscl) {
    x = r_to_x(r);
  } else {
    if (logr) x = log(r);
    else      x = r;
  }

  if (x<xmin) {
    i1 = 0;
    i2 = 1;
  } else if (x>=xmax) {
    i1 = numr - 2;
    i2 = numr - 1;
  } else {
    i1 = (x - xmin)/dx;
    i2 = i1 + 1;
  }

  i1 = std::max<int>(i1, 0);	// Confidence checks
  i2 = std::min<int>(i2, numr-1);

  double A = (xmin + dx*i2 - x)/dx;
  double B = (x - xmin - dx*i1)/dx;

  // Evaluation
  //
  double pot = 0.0, fr = 0.0, ft = 0.0, fp = 0.0;
  for (int l=0; l<=lmax; l+=2) {
    int ll = l/2;

    // integrate positive m values
    for (int m=0; m<=l; m++) {
      double Term1 = T1[ll][m][i1]*A + T1[ll][m][i2]*B;
      double Term2 = T2[ll][m][i1]*A + T2[ll][m][i2]*B;
      double yfac  = Ylm(l, m, cosx)/(2.0*l + 1.0);
      double dfac  = Zlm(l, m, cosx)/(2.0*l + 1.0);

      pot += yfac *      cos(m*phi) * r * (Term1 + Term2);
      fr  += yfac *      cos(m*phi)     * (-Term1*(l+1) + Term2*l);
      fp  += yfac * m * -sin(m*phi) * r * (Term1 + Term2);
      ft  += dfac *      cos(m*phi) * r * (Term1 + Term2);
    }
  }

  pot *= -4.0*M_PI;
  fr  *=  4.0*M_PI;
  ft  *=  4.0*M_PI;
  fp  *=  4.0*M_PI;

  double FR = fr * R/r + ft * z/(r*r);
  double Fz = fr * z/r - ft * R/(r*r);

  return std::tuple<double, double, double, double>(pot, FR, Fz, fp);
}

void DiskEval::write_cache()
{
  std::ofstream cache(cachefile);
  if (cache) {

    // This is a node of simple {key: value} pairs.  More general
    // content can be added as needed.

    YAML::Node node;

    node["model"  ] = model->getID();
    node["mass"   ] = model->getMass();
    if (model->getParams().size())
    node["params" ] = model->getParams();
    node["xmin"   ] = xmin;
    node["xmax"   ] = xmax;
    node["dx"     ] = dx;
    node["ascl"   ] = ascl;
    node["lmax"   ] = lmax;
    node["numr"   ] = numr;
    node["logr"   ] = logr;
    node["xscl"   ] = xscl;

    // Serialize the node
    //
    YAML::Emitter y; y << node;

    // Get the size of the string
    //
    unsigned int hsize = strlen(y.c_str());

    // Write magic #
    //
    cache.write(reinterpret_cast<const char *>(&hmagic),   sizeof(unsigned int));

    // Write YAML string size
    //
    cache.write(reinterpret_cast<const char *>(&hsize),    sizeof(unsigned int));
    
    // Write YAML string
    //
    cache.write(reinterpret_cast<const char *>(y.c_str()), hsize);

    // Cache Term1 and Term2 from the multipole expansion
    //

    // write all m terms (even if not specified)
    //   should record somewhere what mmax is?
    //
    for (int ll=0; ll<=lmax; ll++) {
      for (int mm=0; mm<=lmax; mm++) {
	for (int rr=0; rr<numr; rr++) cache.write((const char *)&T1[ll][mm][rr], sizeof(double));
      }
    }
	
    for (int ll=0; ll<=lmax; ll++) {
      for (int mm=0; mm<=lmax; mm++) {
	for (int rr=0; rr<numr; rr++) cache.write((const char *)&T2[ll][mm][rr], sizeof(double));
      }
    }
    
  } else {
    std::cerr << std::endl
	      << "DiskEval: could not open cache file <"
	      << cachefile << "> for writing" << std::endl;
  }

  std::cerr << std::endl
	    << "DiskEval: wrote cache file <"
	    << cachefile << ">" << std::endl;
}

bool DiskEval::read_cache()
{
  std::ifstream cache(cachefile);
  if (cache) {
    std::string model1;
    double mass1, xmin1, xmax1, dx1, ascl1;
    std::vector<double> params1;
    int lmax1, numr1;
    unsigned char clogr, cxscl;
    bool logr1 = false, xscl1 = false;

    // Attempt to read magic number
    //
    unsigned int tmagic;
    cache.read(reinterpret_cast<char*>(&tmagic), sizeof(unsigned int));

    if (tmagic == hmagic) {
      // YAML size
      //
      unsigned ssize;
      cache.read(reinterpret_cast<char*>(&ssize), sizeof(unsigned int));

      // Make and read char buffer
      //
      auto buf = std::make_unique<char[]>(ssize+1);
      cache.read(buf.get(), ssize);
      buf[ssize] = 0;		// Null terminate

      YAML::Node node;
      
      try {
	node = YAML::Load(buf.get());
      }
      catch (YAML::Exception& error) {
	if (myid)
	  std::cerr << std::endl
		    << "YAML: error parsing <" << buf.get() << "> "
		    << "in " << __FILE__ << ":" << __LINE__ << std::endl
		    << "YAML error: " << error.what() << std::endl;
	throw error;
      }

      // Get parameters
      //
      try {
	model1  = node["model"  ].as<std::string>();
	mass1   = node["mass"   ].as<double>();
	if (node["params"])
	params1 = node["params" ].as<std::vector<double>>();
	xmin1   = node["xmin"   ].as<double>();
	xmax1   = node["xmax"   ].as<double>();
	dx1     = node["dx"     ].as<double>();
	lmax1   = node["lmax"   ].as<int>();
	numr1   = node["numr"   ].as<int>();
	logr1   = node["logr"   ].as<bool>();
	xscl1   = node["xscl"   ].as<bool>();
      }
      catch (YAML::Exception& error) {
	if (myid)
	  std::cerr << "YAML: error reading parameters" << std::endl
		    << "YAML error: " << error.what() << std::endl;
	throw error;
      }

      // Check parameters
      //
      bool okay = true;
   
      if (model1.compare(model->getID()))  {
	okay = false;
	std::cout << std::endl
		  << "DiskEval:read_cache: model ID mismatch <"
		  << model1 << "> != <" << model->getID() << ">"
		  << std::endl;
      }

      // Get parameters
      auto params = model->getParams();

      if (fabs(mass1 - model->getMass()) > 1.0e-18) {
	okay = false;
	std::cout << std::endl
		  << "DiskEval:read_cache: model mass mismatch <"
		  << mass1 << "> != <" << model->getMass() << ">"
		  << std::endl;
      }
    

      // First check parameter size
      if (params1.size() == params.size()) {
	// Now check each parameter
	for (int i=0; i<params.size(); i++) {
	  if (fabs(params1[i] - params[i]) > 1.0e-18) {
	    okay = false;
	    std::cout << std::endl
		      << "DiskEval:read_cache: model parameter ("
		      << i+1 << ") mismatch<" << params1[i] << "> != <"
		      << params[i] << ">" << std::endl;
	  }
	}
      } else {
	    okay = false;
	    std::cout << std::endl
		      << "DiskEval:read_cache: parameter size mismatch<"
		      << params1.size() << "> != <" << params.size() << ">"
		      << std::endl;
      }

      if (fabs(dx1 - dx)     > 1.0e-18) {
	okay = false;
	std::cout << std::endl
		  << "DiskEval:read_cache: dx mismatch <"
		  << dx1 << "> != <" << dx << ">"
		  << std::endl;
      }
    
      if (lmax1 != lmax) {
	okay = false;
	std::cout << std::endl
		  << "DiskEval:read_cache: lmax mismatch <"
		  << lmax1 << "> != <" << lmax << ">"
		  << std::endl;
      }
    
      if (numr1 != numr) {
	okay = false;
	std::cout << std::endl
		  << "DiskEval:read_cache: numr mismatch <"
		  << numr1 << "> != <" << numr << ">"
		  << std::endl;
      }
    
      if ((logr1 and not logr) or (not logr1 and logr)) {
	okay = false;
	std::cout << std::endl
		  << "DiskEval:read_cache: logr mismatch <"
		  << std::boolalpha << logr1 << "> != <"
		  << std::boolalpha << logr << ">" << std::endl;
	
      }
      
      if ((xscl1 and not xscl) or (not xscl1 and xscl)) {
	okay = false;
	std::cout << std::endl
		  << "DiskEval:read_cache: xscl mismatch <"
		  << std::boolalpha << xscl1 << "> != <"
		  << std::boolalpha << xscl << ">" << std::endl;
      }
      
      if (not okay) return false;
    
      // Read multipole expansion terms
      //
      T1.resize(lmax+1);
      T2.resize(lmax+1);

      for (int ll=0; ll<=lmax; ll++) {
	T1[ll].resize(lmax+1);
	T2[ll].resize(lmax+1);
	for (int mm=0; mm<=lmax; mm++) {
	  T1[ll][mm].resize(numr);
	  T2[ll][mm].resize(numr);
	}
      }
      
      for (int ll=0; ll<=lmax; ll++) {
	for (int mm=0; mm<=lmax; mm++) {
	  for (int rr=0; rr<numr; rr++) cache.read((char *)&T1[ll][mm][rr], sizeof(double));
	}
      }
      
      for (int ll=0; ll<=lmax; ll++) {
	for (int mm=0; mm<=lmax; mm++) {
	  for (int rr=0; rr<numr; rr++) cache.read((char *)&T2[ll][mm][rr], sizeof(double));
	}
      }
    }
    
  } else {
    std::cerr << "DiskEval: could not open cache file <"
	      << cachefile << "> for reading" << std::endl
	      << std::endl;
    return false;
  }

  return true;
}

