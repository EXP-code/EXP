#ifndef _SphereEJCOM_H
#define _SphereEJCOM_H

#include <SphereTwoCenter.H>

/** Computes the potential, acceleration and density using the Sturm
    Liouville direct solution fitting to peak potential and COM.
    The mixture function is:
    \f[
    f({\bf r}) = \hbox{erf}\left(C\left|{{\bf r} - {\bf r}_{peak}|\over
    |{\bf r}_{COM} - {\bf r}_{peak}}\right|^\alpha\right)
    \f]
    On evaluation, the acceleration from the grid centered on 
    \f${\bf r}_{COM}\f$ is  multiplied by \f$f({\bf r})\f$ and the
    the acceleration from the grid centered on 
    \f${\bf r}_{peak}\f$ is  multiplied by \f$1-f({\bf r})\f$

    See SphereTwoCenter for other parameters controlling the spherical
    expansion.

    \param cfac is a prefactor $C$ for the erf ramp above
    \param alpha is the exponent in the argument of the erf

*/
class SphereEJCOM : public SphereTwoCenter
{

private:

  double cfac;
  double alpha;

  void initialize();

public:
				// Global parameters
  //! Constructor
  /*! <param> line passes any parameters (rsphSL, rmin and numr) */
  SphereEJCOM(string& line);

  //! Mixture function for the two expansions
  double mixture(Particle& p);

};

#endif


