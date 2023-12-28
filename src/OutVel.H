#ifndef _OutVel_H
#define _OutVel_H

#include <functional>
#include <string>
#include <vector>

#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Attribute.hpp>

#include <interp.H>
#include <Component.H>
#include <EXPException.H>
#include <Coefficients.H>
#include <VelocityBasis.H>

/** Dump coefficients at each interval
    
    @param filename of the coefficient file

    @param name of the component to dump
  
    @param nint is the frequency between file updates 

    @param dof is the number of degrees of freedom. 2=polar, 3 =spherical
    
    @param lmax is the maximum harmonic order

    @param nmax is the maximum radial order

    @param scale is the scale parameter for the density field

    @param rmin is the minimum radius for the density field

    @param rmax is the maximum radius for the density field

    @param ascl is the exponential scale length used for the "expon" density

    @param delta is the cutoff radius used for the "expon" density

    @param model defines the density for the orthogonal functions
*/
class OutVel : public Output
{

private:

  std::string filename, model, outfile;
  Component *tcomp;
  CoefClasses::CoefsPtr coefs;
  int dof;

  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

  using CoefType = Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 3>>;

  //! The basis isntance
  std::shared_ptr<BasisClasses::VelocityBasis> basis;
    
  //! Initialize and write the HDF5 file
  void WriteH5Coefs();

  //! Extend the HDF5 file
  void ExtendH5Coefs();

  //! Write the parameters to the HDF5 file
  void WriteH5Params(HighFive::File& file);

  //! Write the coefficient data to the HDF5 file
  unsigned WriteH5Times(HighFive::Group& snaps, unsigned count);

public:

  //! Constructor
  OutVel(const YAML::Node& conf);

  //! Generate the output
  /*!
    \param nstep is the current time step used to decide whether or not
    to dump
    \param last should be true on final step to force phase space dump
    indepentently of whether or not the frequency criterion is met
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
