#ifndef _EJcom_H
#define _EJcom_H

#include <TwoCenter.H>

/** Computes the potential, acceleration and density using basis
    expansion solutions fitting to peak potential and COM.  This
    implements a TwoCenter class by providing a mixture function.

    The mixture function is:
    \f[
    f({\bf r}) = \hbox{erf}\left(C\left[{|{\bf r} - {\bf r}_{peak}|\over
    |{\bf r}_{COM} - {\bf r}_{peak}|}\right]^\alpha\right)
    \f]
    On evaluation, the acceleration from the grid centered on 
    \f${\bf r}_{COM}\f$ is  multiplied by \f$f({\bf r})\f$ and the
    the acceleration from the grid centered on 
    \f${\bf r}_{peak}\f$ is  multiplied by \f$1-f({\bf r})\f$

    See SphereTwoCenter for other parameters controlling the spherical
    expansion.

    \param cfac is a prefactor \f$C\f$ for the erf ramp above
    \param alpha is the exponent in the argument of the erf (\f$\alpha\f$)

*/
class EJcom : public TwoCenter
{

private:

  //! Prefactor for the erf ramp
  double cfac;

  //! exponent in the argument of the erf
  double alpha;

  //! Initialize the force
  void initialize();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:
				// Global parameters
  /** Constructor
      \param c0 is the instantiating caller (a Component)
      \param line passes any parameters (rsphSL, rmin and numr) 
  */
  EJcom(Component* c0, const YAML::Node& conf);

  //! Mixture function for the two expansions
  double mixture(double* p);

};

#endif
