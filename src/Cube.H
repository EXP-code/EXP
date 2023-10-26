#ifndef _Cube_H
#define _Cube_H

#include <complex>
#include <vector>
#include <string>
#include <set>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor> // For 3d rectangular grids
#include <PotAccel.H>
#include <Coefficients.H>

//! Periodic cube basis
class Cube : public PotAccel
{
private:

  //! Coefficient container
  std::vector<Eigen::Tensor<std::complex<double>, 3>> expccof;

  //! Index extends and particle counts for methods
  int imx, imy, imz, use1, use0;

  //@{
  //! Wave number factors
  double dfac;
  std::complex<double> kfac;
  //@}

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

protected:

  //! Initialize parameters
  void initialize(void);
				// Parameters

  //! Minimum wave numbers in each dimesion
  //@{
  int nmaxx, nmaxy, nmaxz;
  //@}

  //! Maximum wave numbers in each dimesion
  //@{
  int nminx, nminy, nminz;
  //@}

  //@{
  //! Threading
  void * determine_coefficients_thread(void * arg);
  void * determine_acceleration_and_potential_thread(void * arg);
  //@}

  //! Coefficient container instance for writing HDF5
  CoefClasses::CubeCoefs cubeCoefs;

public:

  //! Id string
  std::string id;

  //! Constructor
  Cube(Component* c0, const YAML::Node& conf);

  //! Destructor
  virtual ~Cube();

  //! Compute the coefficients
  void determine_coefficients(void);

  //! Compute the force
  void get_acceleration_and_potential(Component*);

  //! Native binary coefficient dump not implemented
  void dump_coefs(std::ostream& out)
  {
    throw std::runtime_error("Cube: no native binary coefficient dump");
  }

  //! Dump current coefficients into named HDF5 file
  void dump_coefs_h5(const std::string& file);
};


#endif
