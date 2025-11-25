#include <filesystem>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <memory>
#include <limits>
#include <string>

#include <yaml-cpp/yaml.h>	// YAML support

#include <Eigen/Eigen>		// Eigen3
#include <unsupported/Eigen/CXX11/Tensor> // Eigen3 tensors

#include "Progress.H"		// Progress bar

#ifdef HAVE_OMP_H
#include <omp.h>		// For multithreading basis construction
#endif

#include "BiorthCube.H"		// Definition for this class
#include "EXPException.H"	// For GenericError
#include "numerical.H"
#include "libvars.H"

using namespace __EXP__;	// For reference to n-body globals

// Constructor
BiorthCube::BiorthCube(const YAML::Node& conf) : conf(conf)
{
  // Read and assign parameters
  //
  try {
    if (conf["nminx"])       nmin[0] = conf["nminx"].as<int>();
    else                     nmin[0] = std::numeric_limits<int>::max();

    if (conf["nminy"])       nmin[1] = conf["nminy"].as<int>();
    else                     nmin[1] = std::numeric_limits<int>::max();

    if (conf["nminz"])       nmin[2] = conf["nminz"].as<int>();
    else                     nmin[2] = std::numeric_limits<int>::max();

    if (conf["nmaxx"])       nmax[0] = conf["nmaxx"].as<int>();
    else                     nmax[0] = 6;

    if (conf["nmaxy"])       nmax[1] = conf["nmaxy"].as<int>();
    else                     nmax[1] = 6;

    if (conf["nmaxz"])       nmax[2] = conf["nmaxz"].as<int>();
    else                     nmax[2] = 6;

    if (conf["knots"])       knots = conf["knots"].as<int>();
    else                     knots = 24;

    if (conf["verbose"])     verbose = true;
    else                     verbose = false;
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in BiorthCube: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf
			   << std::string(60, '-') << std::endl;
    throw std::runtime_error("BiorthCube: YAML parsing error");
  }

  geometry = "cube";
  forceID  = "BiorthCube";
}


// Get potential basis function
std::complex<double> BiorthCube::pot(Eigen::Vector3d x, Eigen::Vector3i n)
{
  // Skip the constant term
  if (n(0)==0 and n(1)==0 and n(2)==0) return {0};

  double norm = 0.0, sum = 0.0;
  for (int i=0; i<3; i++) {
    norm += n(i)*n(i);
    sum  += x(i)*n(i);
  }

  return std::exp(kfac*sum)/sqrt(nfac*norm);
}

// Get density basis function
std::complex<double> BiorthCube::dens(Eigen::Vector3d x, Eigen::Vector3i n)
{
  // Skip the constant term
  if (n(0)==0 and n(1)==0 and n(2)==0) return {0};

  double norm = 0.0, sum = 0.0;
  for (int i=0; i<3; i++) {
    norm += n(i)*n(i);
    sum  += x(i)*n(i);
  }

  return -std::exp(kfac*sum)*sqrt(nfac*norm);
}

// Get radial force ffrom basis function
Eigen::Vector3cd BiorthCube::force(Eigen::Vector3d x, Eigen::Vector3i n)
{
  double norm = 0.0, sum = 0.0;
  for (int i=0; i<3; i++) {
    norm += n(i)*n(i);
    sum  += x(i)*n(i);
  }

  std::complex<double> fac = std::exp(kfac*sum)/sqrt(nfac*norm);
  Eigen::Vector3cd force;
  for (int i=0; i<3; i++) force(i) = -dfac*n(i)*fac;

  return force;
}

// Get potential
//
std::complex<double>
BiorthCube::get_pot(const BiorthCube::coefType& c, Eigen::Vector3d x)
{
  std::complex<double> potl = 0.0;
    
  // Recursion multipliers and initial values
  Eigen::Vector3cd step, init, curr;
  Eigen::Vector3i i;

  for (int k=0; k<3; k++) {
    step(k) = std::exp(kfac*x(k));
    init(k) = std::exp(-kfac*(x(k)*nmax(k)));
  }
    
  curr(0) = init(0);
  for (i(0)=0; i(0)<=2*nmax(0); i(0)++, curr(0)*=step(0)) {
    curr(1) = init(1);
    for (i(1)=0; i(1)<=2*nmax(1); i(1)++, curr(1)*=step(1)) {
      curr(2) = init(2);
      for (i(2)=0; i(2)<=2*nmax(2); i(2)++, curr(2)*=step(2)) {
	
	auto fac = curr(0)*curr(1)*curr(2)*c(i(0), i(1), i(2));
	  
	// Compute wavenumber; the coefficients are stored as:
	// -nmax,-nmax+1,...,0,...,nmax-1,nmax
	//
	Eigen::Vector3i ii = i - nmax;

	// No contribution to acceleration and potential ("swindle")
	// for zero wavenumber
	if (ii(0)==0 && ii(1)==0 && ii(2)==0) continue;
	  
	// Limit to minimum wave number
	if (abs(ii(0)) > nmin(0) ||
	    abs(ii(1)) > nmin(1) ||
	    abs(ii(2)) > nmin(2)  ) continue;
	  
	// Normalization
	double norm = 1.0/sqrt(nfac*ii.dot(ii));

	potl += fac*norm;
      }
    }
  }

  return potl;
}

// Density evaluation
//
std::complex<double>
BiorthCube::get_dens(const BiorthCube::coefType& c, Eigen::Vector3d x)
{
  const double dfac = 2.0*M_PI;
  const std::complex<double> kfac = std::complex<double>(0.0, dfac);
  std::complex<double> dens = 0.0;
    
  // Recursion multipliers and initial values
  Eigen::Vector3cd step, init, curr;
  Eigen::Vector3i i;

  for (int k=0; k<3; k++) {
    step(k) = std::exp(kfac*x(k));
    init(k) = std::exp(-kfac*(x(k)*nmax(k)));
  }
    
  curr(0) = init(0);
  for (i(0)=0; i(0)<=2*nmax(0); i(0)++, curr(0)*=step(0)) {
    curr(1) = init(1);
    for (i(1)=0; i(1)<=2*nmax(1); i(1)++, curr(1)*=step(1)) {
      curr(2) = init(2);
      for (i(2)=0; i(2)<=2*nmax(2); i(2)++, curr(2)*=step(2)) {
	
	auto fac = curr(0)*curr(1)*curr(2)*c(i(0), i(1), i(2));
	  
	// Compute wavenumber; the coefficients are stored as:
	// -nmax,-nmax+1,...,0,...,nmax-1,nmax
	//
	Eigen::Vector3i ii = i - nmax;

	// No contribution to acceleration and potential ("swindle")
	// for zero wavenumber
	if (ii(0)==0 && ii(1)==0 && ii(2)==0) continue;
	  
	// Limit to minimum wave number
	// Limit to minimum wave number
	if (abs(ii(0)) > nmin(0) ||
	    abs(ii(1)) > nmin(1) ||
	    abs(ii(2)) > nmin(2)  ) continue;
	  
	// Normalization
	double norm = -sqrt(nfac*ii.dot(ii));

	dens += fac*norm;
      }
    }
  }

  return dens;
}


// Force evaluation
Eigen::Vector3cd
BiorthCube::get_force(const BiorthCube::coefType& c, Eigen::Vector3d x)
{
  const double dfac = 2.0*M_PI;
  const std::complex<double> kfac = std::complex<double>(0.0, dfac);
  Eigen::Vector3cd force{0.0, 0.0, 0.0};
  
    
  // Recursion multipliers and initial values
  Eigen::Vector3cd step, init, curr;
  Eigen::Vector3i i;

  for (int k=0; k<3; k++) {
    step(k)  = std::exp(kfac*x(k));
    init(k) = std::exp(-kfac*(x(k)*nmax(k)));
  }
    
  curr(0) = init(0);
  for (i(0)=0; i(0)<=2*nmax(0); i(0)++, curr(0)*=step(0)) {
    curr(1) = init(1);
    for (i(1)=0; i(1)<=2*nmax(1); i(1)++, curr(1)*=step(1)) {
      curr(2) = init(2);
      for (i(2)=0; i(2)<=2*nmax(2); i(2)++, curr(2)*=step(2)) {
	
	auto fac = curr(0)*curr(1)*curr(2)*c(i(0), i(1), i(2));
	  
	// Compute wavenumber; the coefficients are stored as:
	// -nmax,-nmax+1,...,0,...,nmax-1,nmax
	//
	Eigen::Vector3i ii = i - nmax;

	// No contribution to acceleration and potential ("swindle")
	// for zero wavenumber
	if (ii(0)==0 && ii(1)==0 && ii(2)==0) continue;
	  
	// Limit to minimum wave number
	if (abs(ii(0)) > nmin(0) ||
	    abs(ii(1)) > nmin(1) ||
	    abs(ii(2)) > nmin(2)  ) continue;
	  
	// Normalization
	double norm = 1.0/sqrt(nfac*ii.dot(ii));

	for (int k=0; k<3; k++)
	  force(k) -= std::complex<double>(0.0, dfac*i(k))*fac*norm;
      }
    }
  }

  return force;
}

Eigen::MatrixXcd BiorthCube::orthoCheck()
{
  Eigen::Vector3i d {2*nmax(0)+1, 2*nmax(1)+1, 2*nmax(2)+1};
  int dim = d(0)*d(1)*d(2);
  Eigen::MatrixXcd ortho(dim, dim);
  ortho.setZero();

  double dx = 1.0/knots, dy = 1.0/knots, dz = 1.0/knots;
  double vol = dx*dy*dz;

  std::shared_ptr<progress::progress_display> progress;
  if (myid==0 and verbose)
    progress = std::make_shared<progress::progress_display>(knots*knots*knots);

  auto index = [&](int indx)
  {
    Eigen::Vector3i n;

    // Flattened index is: n0*d1*d2 + n1*d2 + n2

    n(0) = indx/d(1)/d(2);	  // indx/(d1*d2)
    n(1) = indx/d(2) - n(0)*d(1); // (indx - n0*d1*d2)/d2
    n(2) = indx - (n(0)*d(1) + n(1))*d(2); // indx - n0*d1*d2 - n1*d2
    n -= nmax;

    return n;
  };


  // Centered rectangle quadrature for Fourier basis
  //
  Eigen::Vector3d x;		// Position vector
  for (int i=0; i<knots; i++) {
    x(0) = dx*i;

    for (int j=0; j<knots; j++) {
      x(1) = dy*j;

      for (int k=0; k<knots; k++) {
	x(2) = dz*k;

#pragma omp parallel for
	for (int indx=0; indx<dim*dim; indx++) {
	  int indx1 = indx/dim;
	  int indx2 = indx - indx1*dim;

	  ortho(indx1, indx2) += vol *
	    pot(x, index(indx1)) * std::conj(-dens(x, index(indx2)));
	}
	// Progress bar
	if (progress) ++(*progress);
      }
    }
  }
  
  // Artifically insert the undefined constant term
  //
  int zindx = nmax(0)*(2*nmax(1)+1)*(2*nmax(2)+1) + nmax(1)*(2*nmax(2)+1) + nmax(2);
  auto nn = index(zindx);
  ortho(zindx, zindx) = 1.0;

  // Internal diagnostic check
  //
  if (verbose) {
    double maxdiag = 0.0;
    double maxoffd = 0.0;
    std::map<double, std::pair<int, int>> worst;
    for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
	double bad;
	if (i==j) {
	  bad = std::abs(ortho(i, j) - 1.0);
	  maxdiag = std::max<double>(maxdiag, bad);
	} else {
	  bad = std::abs(ortho(i, j));
	  maxoffd = std::max<double>(maxoffd, bad);
	}
	worst[bad] = {i, j};
      }
    }
    
    std::cout << "Orthogonality check:"
	      << " diag=" << maxdiag
	      << " offd=" << maxoffd << std::endl << std::endl;
    const int shame = 10;
    std::cout << "Top " << shame << " offenders:" << std::endl;
    auto rit = worst.rbegin();
    for (int i=0; i<shame; i++, rit++) {
      auto n1 = index(rit->second.first);
      auto n2 = index(rit->second.second);
      std::ostringstream s1, s2;
      s1 << '[' << n1(0) << ',' << n1(1) << ',' << n1(2) << ']';
      s2 << '[' << n2(0) << ',' << n2(1) << ',' << n2(2) << ']';
      std::cout << "** " << std::setw(18) << rit->first
		<< std::setw(16) << s1.str()
		<< std::setw(16) << s2.str()
		<< std::setw( 8) << rit->second.first
		<< std::setw( 8) << rit->second.second
		<< std::setw( 8) << dim
		<< std::endl;
    }
  }

  return ortho;
}
