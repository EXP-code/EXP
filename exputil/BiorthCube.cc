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

#ifdef HAVE_OMP_H
#include <omp.h>		// For multithreading basis construction
#endif

#include <BiorthCube.H>		// Definition for this class
#include <EXPException.H>	// For GenericError
#include <numerical.H>
#include <libvars.H>

using namespace __EXP__;	// For reference to n-body globals

// Constructor
BiorthCube::BiorthCube(const YAML::Node& conf) : conf(conf)
{
  // Read and assign parameters
  //
  try {
    if (conf["nminx"])       nminx = conf["nminx"].as<int>();
    else                     nminx = std::numeric_limits<int>::max();

    if (conf["nminy"])       nminy = conf["nminy"].as<int>();
    else                     nminy = std::numeric_limits<int>::max();

    if (conf["nminz"])       nminz = conf["nminz"].as<int>();
    else                     nminz = std::numeric_limits<int>::max();

    if (conf["nmaxx"])       nmaxx = conf["nmaxx"].as<int>();
    else                     nmaxx = 6;

    if (conf["nmaxy"])       nmaxy = conf["nmaxy"].as<int>();
    else                     nmaxy = 6;

    if (conf["nmaxz"])       nmaxz = conf["nmaxz"].as<int>();
    else                     nmaxz = 6;

    if (conf["knots"])       knots = conf["knots"].as<int>();
    else                     knots = 24;
  }
  catch (YAML::Exception & error) {
    if (myid==0) std::cout << "Error parsing parameters in BiorthCube: "
			   << error.what() << std::endl
			   << std::string(60, '-') << std::endl
			   << "Config node"        << std::endl
			   << std::string(60, '-') << std::endl
			   << conf
			   << std::string(60, '-') << std::endl;
    MPI_Finalize();
    exit(-1);
  }

  geometry = "cube";
  forceID  = "BiorthCube";
}


// Get potential basis function
std::complex<double> BiorthCube::pot(Eigen::Vector3d x, Eigen::Vector3i n)
{
  double norm = 0.0, sum = 0.0;
  for (int i=0; i<3; i++) {
    norm += n(i)*n(i);
    sum  += x(i)*n(i);
  }

  std::complex<double> rad{0.0, 2.0*M_PI*sum};
  return std::exp(rad)/sqrt(norm);
}

// Get density basis function
std::complex<double> BiorthCube::dens(Eigen::Vector3d x, Eigen::Vector3i n)
{
  double norm = 0.0, sum = 0.0;
  for (int i=0; i<3; i++) {
    norm += n(i)*n(i);
    sum  += x(i)*n(i);
  }

  std::complex<double> rad{0.0, 2.0*M_PI*sum};
  return -std::exp(rad)*sqrt(norm);
}

// Get radial force ffrom basis function
Eigen::Vector3cd BiorthCube::force(Eigen::Vector3d x, Eigen::Vector3i n)
{
  double norm = 0.0, sum = 0.0;
  for (int i=0; i<3; i++) {
    norm += n(i)*n(i);
    sum  += x(i)*n(i);
  }

  std::complex<double> rad{0.0, 2.0*M_PI*sum};
  std::complex<double> fac = std::exp(rad)/sqrt(norm);
  Eigen::Vector3cd force;
  for (int i=0; i<3; i++) force(i) = -2.0*M_PI*n(i)*fac;

  return force;
}

// Get potential
//
std::complex<double>
BiorthCube::get_pot(const BiorthCube::coefType& c, Eigen::Vector3d x)
{
  const double dfac = 2.0*M_PI;
  const std::complex<double> kfac = std::complex<double>(0.0, dfac);
  std::complex<double> potl = 0.0;
    
  // Recursion multipliers and initial values
  Eigen::Vector3cd step, curr;
  for (int i=0; i<3; i++) {
    step(i) = std::exp(kfac*x(i));
    curr(i) = std::exp(-kfac*static_cast<double>(nmaxx)*x(i));
  }
    
  for (int ix=0; ix<=2*nmaxx; ix++, curr(0)*=step(0)) {
    for (int iy=0; iy<=2*nmaxy; iy++, curr(1)*=step(1)) {
      for (int iz=0; iz<=2*nmaxz; iz++, curr(2)*=step(2)) {
	
	auto fac = curr(0)*curr(1)*curr(2)*c(ix, iy, iz);
	  
	// Compute wavenumber; recall that the coefficients are
	// stored as follows: -nmax,-nmax+1,...,0,...,nmax-1,nmax
	//
	int ii = ix - nmaxx;
	int jj = iy - nmaxy;
	int kk = iz - nmaxz;
	  

	// No contribution to acceleration and potential ("swindle")
	// for zero wavenumber
	if (ii==0 && jj==0 && kk==0) continue;
	  
	// Limit to minimum wave number
	if (abs(ii)<nminx || abs(jj)<nminy || abs(kk)<nminz) continue;
	  
	// Normalization
	double norm = 1.0/sqrt(M_PI*(ii*ii + jj*jj + kk*kk));;

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
  Eigen::Vector3cd step, curr;
  for (int i=0; i<3; i++) {
    step(i)  = std::exp(kfac*x(i));
    curr(i) = std::exp(-kfac*static_cast<double>(nmaxx)*x(i));
  }
    
  for (int ix=0; ix<=2*nmaxx; ix++, curr(0)*=step(0)) {
    for (int iy=0; iy<=2*nmaxy; iy++, curr(1)*=step(1)) {
      for (int iz=0; iz<=2*nmaxz; iz++, curr(2)*=step(2)) {
	
	auto fac = curr(0)*curr(1)*curr(2)*c(ix, iy, iz);
	  
	// Compute wavenumber; recall that the coefficients are
	// stored as follows: -nmax,-nmax+1,...,0,...,nmax-1,nmax
	//
	int ii = ix - nmaxx;
	int jj = iy - nmaxy;
	int kk = iz - nmaxz;
	  

	// No contribution to acceleration and potential ("swindle")
	// for zero wavenumber
	if (ii==0 && jj==0 && kk==0) continue;
	  
	// Limit to minimum wave number
	if (abs(ii)<nminx || abs(jj)<nminy || abs(kk)<nminz) continue;
	  
	// Normalization
	double norm = -sqrt(M_PI*(ii*ii + jj*jj + kk*kk));;

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
  Eigen::Vector3cd step, curr;
  for (int i=0; i<3; i++) {
    step(i)  = std::exp(kfac*x(i));
    curr(i) = std::exp(-kfac*static_cast<double>(nmaxx)*x(i));
  }
    
  for (int ix=0; ix<=2*nmaxx; ix++, curr(0)*=step(0)) {
    for (int iy=0; iy<=2*nmaxy; iy++, curr(1)*=step(1)) {
      for (int iz=0; iz<=2*nmaxz; iz++, curr(2)*=step(2)) {
	
	auto fac = curr(0)*curr(1)*curr(2)*c(ix, iy, iz);
	  
	// Compute wavenumber; recall that the coefficients are
	// stored as follows: -nmax,-nmax+1,...,0,...,nmax-1,nmax
	//
	int ii = ix - nmaxx;
	int jj = iy - nmaxy;
	int kk = iz - nmaxz;
	  

	// No contribution to acceleration and potential ("swindle")
	// for zero wavenumber
	if (ii==0 && jj==0 && kk==0) continue;
	  
	// Limit to minimum wave number
	if (abs(ii)<nminx || abs(jj)<nminy || abs(kk)<nminz) continue;
	  
	// Normalization
	double norm = 1.0/sqrt(M_PI*(ii*ii + jj*jj + kk*kk));;

	force(0) -= std::complex<double>(0.0, dfac*ii)*fac*norm;
	force(1) -= std::complex<double>(0.0, dfac*jj)*fac*norm;
	force(2) -= std::complex<double>(0.0, dfac*kk)*fac*norm;
      }
    }
  }

  return force;
}

Eigen::MatrixXcd BiorthCube::orthoCheck()
{
  Eigen::Vector3i d{2*nmaxx+1, 2*nmaxy+1, 2*nmaxz+1};
  int dim = d(0)*d(1)*d(2);
  Eigen::MatrixXcd ortho(dim, dim);
  ortho.setZero();

  double dx = 1.0/knots, dy = 1.0/knots, dz = 1.0/knots;
  double vol = dx*dy*dz;

  Eigen::Vector3d x;
  Eigen::Vector3i n1, n2;

  for (int i=0; i<knots; i++) {
    x(0) = (0.5 + i)*dx;

    for (int j=0; j<knots; j++) {
      x(1) = (0.5 + j)*dy;

      for (int k=0; i<knots; k++) {
	x(2) = (0.5 + k)*dz;

	for (int indx1=0; indx1<dim; indx1++) {

	  n1(0) = indx1/d(1)/d(2);
	  n1(1) = indx1/d(2) - n1(0)*d(1);
	  n1(2) = indx1 - (n1(0)*d(1) + n1(1))*d(2);

	  for (int indx2=0; indx2<dim; indx2++) {

	    n2(0) = indx2/d(1)/d(2);
	    n2(1) = indx2/d(2) - n2(0)*d(1);
	    n2(2) = indx2 - (n2(0)*d(1) + n2(1))*d(2);

	    ortho(indx1, indx2) += vol *
	      pot(x, n1) * std::conj(dens(x, n2));
	      
	  }
	}
      }
    }
  }
  
  return ortho;
}
