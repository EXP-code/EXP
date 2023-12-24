#ifndef _OutVel_H
#define _OutVel_H

#include <functional>
#include <string>
#include <vector>

#include <Eigen/Eigen>

#include <interp.H>
#include <Component.H>
#include <EXPException.H>
#include <OrthoFunction.H>

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

  std::string filename, model;
  Component *tcomp;
  double scale, rmin, rmax, ascl, delta;
  int lmax, nmax, dof;

  std::shared_ptr<Linear1d> interp;
  std::function<double(double)> densfunc;

  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

  //! Coefficient tables
  std::vector<Eigen::MatrixXcd> coefs;

  //! The orthgonal function instance
  std::shared_ptr<OrthoFunction> ortho;

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
