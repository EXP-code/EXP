#include <iostream>
#include <iomanip>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/progress.hpp>	// Progress bar

#include <DiskEval.H>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

const std::string DiskEval::cachefile = ".DiskEval.cache";

DiskEval::DiskEval
(EmpCylSL::AxiDiskPtr model, double rmin, double rmax, double ascl,
 int lmax, int numr, int nint, bool use_progress) :
  model(model), ascl(ascl), lmax(lmax), numr(numr)
{
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

  if (read_cache()) return;

  // First: compute the disk density components
  //
  rho.resize(lmax+1);
  for (auto & v : rho) v.resize(numr, 0.0);
  
  // Gauss-Legendre
  //
  boost::shared_ptr<LegeQuad> lq = boost::make_shared<LegeQuad>(nint);

  int numthrd = 1, tid = 0;
#ifdef HAVE_OPENMP
#pragma omp parallel
  {
    numthrd = omp_get_num_threads();
    int myid = omp_get_thread_num();
    if (myid==0)
      std::cout << "Number of threads=" << numthrd << std::endl;
  }
#endif

  boost::shared_ptr<boost::progress_display> progress;
  if (use_progress) {
    std::cout << std::endl << "Begin: exact force evaluation"
	      << std::endl << "-----------------------------"
	      << std::endl;
      
    std::cout << std::endl << "Multipole density coefficients"
	      << std::endl;
    progress = boost::make_shared<boost::progress_display>(numr/numthrd);
  }
  
  // Compute \rho_{lm} on radial grid
  //
#pragma omp parallel for
  for (int i=0; i<numr; i++) {
#ifdef HAVE_OPENMP
    tid = omp_get_thread_num();
#endif

    double x = xmin + dx*i, r;

    if (xscl) {
      r = x_to_r(x);
    } else {
      if (logr) r = exp(x);
      else      r = x;
    }

    for (int n=1; n<=nint; n++) {
      double cosx = lq->knot(n); // Assume rho(R, z) = rho(R, -z);
      double R = sqrt(1.0 - cosx*cosx) * r;
      double z = cosx * r;

      double dens = (*model)(R, z);

      double fac = 2.0*M_PI * 2.0*lq->weight(n);

      // Even terms only--------+
      //                        |
      //                        v
      for (int l=0; l<=lmax; l+=2) {
	int ll = l/2;
	rho[ll][i] += Ylm(l, 0, cosx) * dens * fac;
      }
    }

    // Progress bar
    //
    if (progress and tid==0) ++(*progress);
  }

  // Compute Term1 and Term2 by quadrature
  //
  T1.resize(lmax+1);
  for (auto & v : T1) v.resize(numr);

  T2.resize(lmax+1);
  for (auto & v : T2) v.resize(numr);
  

  int numrads = 0;
  for (int l=0; l<=lmax; l+=2) numrads++;
  numrads *= numr;
  numrads /= numthrd;

  if (use_progress) {
    std::cout << std::endl << "Quadrature loop multipole expansion"
	      << std::endl;
    progress = boost::make_shared<boost::progress_display>(numrads);
  }
  
  // l loop
  //
  for (int l=0; l<=lmax; l+=2) {
    int ll = l/2;

    // Outer r loop
    //
#pragma omp parallel for
    for (int i=0; i<numr; i++) {
#ifdef HAVE_OPENMP
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
	    sum1 += 0.5*(rho[ll][j-1]*pow(rl/rr, l+2)*dr_to_dx(xl) +
			 rho[ll][j  ]*pow(rp/rr, l+2)*dr_to_dx(xp)) * dx;
	  } else {
	    if (logr) {
	      sum1 += 0.5*(rho[ll][j-1]*pow(rl/rr, l+3) +
			   rho[ll][j  ]*pow(rp/rr, l+3)) * rr * dx;
	    } else {
	      sum1 += 0.5*(rho[ll][j-1]*pow(rl/rr, l+2) +
			   rho[ll][j  ]*pow(rp/rr, l+2)) * dx;
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
	  sum2 += 0.5*(rho[ll][j  ]*pow(rl/rr, 1-l) * dr_to_dx(xl) +
		       rho[ll][j+1]*pow(rp/rr, 1-l) * dr_to_dx(xp)) * dx;
	} else {
	  if (logr) {
	    sum2 += 0.5*(rho[ll][j  ]*pow(rl/rr, 2-l) +
			 rho[ll][j+1]*pow(rp/rr, 2-l)) * rr * dx;
	  } else {
	    sum2 += 0.5*(rho[ll][j  ]*pow(rl/rr, 1-l) +
			 rho[ll][j+1]*pow(rp/rr, 1-l)) * dx;
	  }
	}
      }

      // Save integral values
      //
      T1[ll][i] = sum1;
      T2[ll][i] = sum2;

      // Progress bar
      //
      if (progress and tid==0) ++(*progress);
    }
    // END: outer loop over r

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
	test << std::setw(18) << rho[ll][i];
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
	test << std::setw(18) << T1[ll][i];
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
	test << std::setw(18) << T2[ll][i];
      }
      test << std::endl;
    }
  }

  write_cache();
}


std::tuple<double, double, double> DiskEval::operator()(double R, double z)
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

  i1 = std::max<int>(i1, 0);	// Sanity checks
  i2 = std::min<int>(i2, numr-1);

  double A = (xmin + dx*i2 - x)/dx;
  double B = (x - xmin - dx*i1)/dx;

  // Evaluation
  //
  double pot = 0.0, fr = 0.0, ft = 0.0;
  for (int l=0; l<=lmax; l+=2) {
    int ll = l/2;
    double Term1 = T1[ll][i1]*A + T1[ll][i2]*B;
    double Term2 = T2[ll][i1]*A + T2[ll][i2]*B;
    double yfac  = Ylm(l, 0, cosx)/(2.0*l + 1.0);
    double dfac  = Zlm(l, 0, cosx)/(2.0*l + 1.0);

    pot += yfac * r * (Term1 + Term2);
    fr  += yfac * (-Term1*(l+1) + Term2*l);
    ft  += dfac * r * (Term1 + Term2);
  }

  pot *= -4.0*M_PI;
  fr  *=  4.0*M_PI;
  ft  *=  4.0*M_PI;

  double FR = fr * R/r + ft * z/(r*r);
  double Fz = fr * z/r - ft * R/(r*r);

  return std::tuple<double, double, double>(pot, FR, Fz);
}

void DiskEval::write_cache()
{
  std::ofstream cache(cachefile);
  if (cache) {
    char str[32];
    std::fill(str, str+32, 0);
    std::string id = model->getID();
    int idlen = id.size();
    strncpy(str, id.c_str(), idlen);
    unsigned char clogr = 0, cxscl = 0;
    if (logr) clogr = 1;
    if (xscl) cxscl = 1;

    // Parameters
    //
    cache.write((const char *)str,    sizeof(char)*32);
    cache.write((const char *)&xmin,  sizeof(double));
    cache.write((const char *)&xmax,  sizeof(double));
    cache.write((const char *)&dx,    sizeof(double));
    cache.write((const char *)&ascl,  sizeof(double));
    cache.write((const char *)&lmax,  sizeof(int   ));
    cache.write((const char *)&numr,  sizeof(int   ));
    cache.write((const char *)&clogr, sizeof(unsigned char));
    cache.write((const char *)&cxscl, sizeof(unsigned char));

    // Cache Term1 and Term2 from the multipole expansion
    //
    for (auto v : T1) {
      cache.write((const char *)v.data(), sizeof(double)*numr);
    }
    for (auto v : T2) {
      cache.write((const char *)v.data(), sizeof(double)*numr);
    }

  } else {
    std::cerr << "DiskEval: could not open cache file <"
	      << cachefile << "> for writing" << std::endl;
  }

  std::cerr << "DiskEval: wrote cache file <"
	    << cachefile << ">" << std::endl;
}

bool DiskEval::read_cache()
{
  std::ifstream cache(cachefile);
  if (cache) {
    char str[32];
    double xmin1, xmax1, dx1, ascl1;
    int lmax1, numr1;
    unsigned char clogr, cxscl;
    bool logr1 = false, xscl1 = false;

    // Parameters
    //
    cache.read((char *)str,    sizeof(char)*32);
    cache.read((char *)&xmin1, sizeof(double));
    cache.read((char *)&xmax1, sizeof(double));
    cache.read((char *)&dx1,   sizeof(double));
    cache.read((char *)&ascl1, sizeof(double));
    cache.read((char *)&lmax1, sizeof(int   ));
    cache.read((char *)&numr1, sizeof(int   ));
    cache.read((char *)&clogr, sizeof(unsigned char));
    cache.read((char *)&cxscl, sizeof(unsigned char));

    if (clogr) logr1 = true;
    if (cxscl) xscl1 = true;

    std::string ID1(str);

    // Check parameters
    bool okay = true;
   
    if (ID1.compare(model->getID()))  {
      okay = false;
      std::cout << "DiskEval:read_cache: model ID mismatch <"
		<< ID1 << "> != <" << model->getID() << ">"
		<< std::endl;
    }
    if (fabs(xmin1 - xmin) > 1.0e-18) {
      okay = false;
      std::cout << "DiskEval:read_cache: xmin mismatch <"
		<< xmin1 << "> != <" << xmin << ">"
		<< std::endl;
    }

    if (fabs(xmax1 - xmax) > 1.0e-18) {
      okay = false;
      std::cout << "DiskEval:read_cache: xmax mismatch <"
		<< xmax1 << "> != <" << xmax << ">"
		<< std::endl;
    }
    
    if (fabs(dx1 - dx)     > 1.0e-18) {
      okay = false;
      std::cout << "DiskEval:read_cache: dx mismatch <"
		<< dx1 << "> != <" << dx << ">"
		<< std::endl;
    }
    
    if (lmax1 != lmax) {
      okay = false;
      std::cout << "DiskEval:read_cache: lmax mismatch <"
		<< lmax1 << "> != <" << lmax << ">"
		<< std::endl;
    }
    
    if (numr1 != numr) {
      okay = false;
      std::cout << "DiskEval:read_cache: numr mismatch <"
		<< numr1 << "> != <" << numr << ">"
		<< std::endl;
    }
    
    if ((logr1 and not logr) or (not logr1 and logr)) {
      okay = false;
      std::cout << "DiskEval:read_cache: logr mismatch <"
		<< std::boolalpha << logr1 << "> != <"
		<< std::boolalpha << logr << ">" << std::endl;
      
    }
    
    if ((xscl1 and not xscl) or (not xscl1 and xscl)) {
      okay = false;
      std::cout << "DiskEval:read_cache: xscl mismatch <"
		<< std::boolalpha << xscl1 << "> != <"
		<< std::boolalpha << xscl << ">" << std::endl;
    }
    
    if (not okay) return false;
    
    // Read multipole expansion terms
    //
    T1.resize(lmax+1);
    T2.resize(lmax+1);
    
    for (auto & v : T1) {
      v.resize(numr);
      cache.read((char *)v.data(), sizeof(double)*numr);
    }
    
    for (auto & v : T2) {
      v.resize(numr);
      cache.read((char *)v.data(), sizeof(double)*numr);
    }
    
  } else {
    std::cerr << "DiskEval: could not open cache file <"
	      << cachefile << "> for reading" << std::endl
	      << std::endl;
    return false;
  }

  return true;
}
