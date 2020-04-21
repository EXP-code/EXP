#include <iostream>
#include <iomanip>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/progress.hpp>	// Progress bar

#include <DiskEval.H>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

DiskEval::DiskEval
(EmpCylSL::AxiDiskPtr model, double rmin, double rmax, double ascl,
 int lmax, int numr, int nint, bool use_progress) :
  model(model), ascl(ascl), lmax(lmax), numr(numr)
{
  if (ascl>0.0) xscl = true;
  else          xscl = false;

  // Assign grid parameters
  //
  if (xscl) {
    xmin = r_to_x(rmin);
    xmax = r_to_x(rmax);
  } else {
    if (rmin > 0.0) {
      xmin = log(rmin);
      xmax = log(rmax);
      logr = true;
    } else {
      xmin = rmin;
      xmax = rmax;
      logr = false;
    }
  }

  dx = (xmax - xmin)/(numr-1);

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

      double xx = xmin + dx*(i-1), rr;

      if (xscl) {
	rr = x_to_r(xx);
      } else {
	if (logr) rr = exp(xx);
	else      rr = xx;
      }

      // Integral over 0 to r
      //
      double sum1 = 0.0, sum2 = 0.0;

      for (int j=1; j<=i; j++) {

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

