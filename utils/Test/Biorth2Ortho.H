#pragma once

#include <cmath>
#include <cstdlib>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>		// For potential error handling
#include <functional>
#include <type_traits>		// For std::enable_if and type traits

#include <Eigen/Eigen>		// Eigen library


// HighFive HDF5 C++ wrapper
#include <highfive/highfive.hpp>
#include <highfive/eigen.hpp>

// EXP framework
#include <Progress.H>
#include <localmpi.H>
#include <biorth.H>

/** Class to convert a biorthogonal basis to an orthonormal basis

    Beginning with a biorthogonal basis set, one of the dual functions
    is used as the support for constrruction an orthonormal basis set
    using Gram-Schmidt orthogonalization.  Transformations between the
    orthonormal and biorthogonal bases are computed and can be used to
    convert coefficients or functions to and from the dual space.
 */
class Biorth2Ortho
{
public:

  using Real = long double;
  using MatrixXld = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
  using VectorXld = Eigen::Vector<Real, Eigen::Dynamic>;

protected:

  std::shared_ptr<AxiSymBiorth> biorth;

  int nmax;			// Radial basis order
  int lmax;			// Angular basis order
  int ngrid;			// Number of radial grid points

  double rmin;			// Minimum expansion radius
  double rmax;			// Maximum expansion radius
  double scl;			// Scale for Laguerre basis

  double xmin;			// Minimum grid value
  double xmax;			// Maximum grid value
  double dx;			// Grid spacing

  bool prog_bar = true;		// Use progress bar

  //! Basis transformation matrices
  std::vector<MatrixXld> Trho, Tphi;

  //! Basis grids
  std::vector<MatrixXld> Wtbl;

  //@{
  //! Algorithm switches
  bool laguerre = false;	// Use Laguerre polynomials as initial basis (for testing)
  bool alternate = false;	// Alternate density and potential functions (for testing; not recommended)
  bool sum = false;		// Use sum of density and potential as initial basis (for testing; not recommended)
  bool classic = false;		// Use classical Gram-Schmidt (else modified Gram-Schmidt; not recommended)
  bool density = true;		// Use density functions (else potential; recommended)
  bool weight = false;	        // Use density weighting in inner product (not recommended)

				// Biorthogonal function type for basis
  ScalarType stype = ScalarType::density;

  int myid = 0;			// MPI rank
  int numprocs = 1;		// Number of MPI processes
  //@}

  class HDF5Exception : public std::runtime_error
  {
  public:
    HDF5Exception(const std::string message) : std::runtime_error(message) {}
  };

  class MismatchException : public std::runtime_error
  {
  public:
    MismatchException(const std::string message) : std::runtime_error(message) {}
  };

public:

  //! Constructor
  Biorth2Ortho(std::shared_ptr<AxiSymBiorth> b, int Lmax, int Nmax, int Ngrid,
	       double Rmin, double Rmax, double scale, bool weight=false);

  //! Destructor
  ~Biorth2Ortho(void);

  //! Convert biorthogonal to orthonormal basis
  void generate();

  //! Turn off progress bar
  void noBar(void) { prog_bar = false; }

  //! Set Laguerre basis
  void useLaguerre(bool tf=true) { laguerre = tf; }

  //! Set alternate density/potential basis
  void useAlternate(bool tf=true) { alternate = tf; }

  //! Set sum basis
  void useSum(bool tf=true) { sum = tf; }

  //! Set classical Gram-Schmidt
  void useClassic(bool tf=true) { classic = tf; }

  //! Basis generation type
  void useDensity(bool tf=true) {
    stype = tf ? ScalarType::density : ScalarType::potential;
  }


  //@{

  //! Access to transformation matrices. Transforms from the
  //! orthonormal basis and the biorthogonal basis.  These are one of
  //! tha main products of this class.
  const std::vector<MatrixXld>& getTrho(void) const { return Trho; }
  const std::vector<MatrixXld>& getTphi(void) const { return Tphi; }
  //@}

  //! Access to new basis functions on grid
  const std::vector<MatrixXld>& getWtbl(void) const { return Wtbl; }

  //! Interpolate on the orthornomal basis functions
  const double basis(int l, int n, double r);

  //! Write output file with transformation matrices and diagnostics
  void output(const std::string& prefix);

  //! Write an HDF5 cache file
  void writeH5(const std::string& filename);

  //! Read an HDF5 cache file
  void readH5(const std::string& filename, bool check=true);

};
