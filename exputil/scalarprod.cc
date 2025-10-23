/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Determine inner product for a given function
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 02/15/94
 *
 ***************************************************************************/


#include <functional>
#include <cstdlib>
#include <cmath>

#include "gaussQ.H"
#include "biorth2d.H"
#include "OrthoPoly.H"


Eigen::VectorXd
scalar_prod(ScalarType type, double rmin, double rmax, int l, int m,
	    AxiSymBiorth& s, std::function<double(double, int, int)> func, 
	    int numc, int numg)
{

  if ( !(s.get_dof()==2 || s.get_dof()==3) ) {
    cerr << "scalar_prod: dof=" << s.get_dof() << ",  must be 2 or 3!\n";
    exit(-1);
  }

  LegeQuad qe(numg);

  Eigen::VectorXd coef(numc);
  coef.setZero();

  double r, del = rmax - rmin;

  for (int i=0; i<numg; i++) {
    r = rmin + del*qe.knot(i);
    for (int n=0; n<numc; n++) {
      switch (type) {

      case density:
	if (s.get_dof()==2)
	  coef[n] += del * qe.weight(i) * pow(r, (int)(s.get_dof()-1)) *
	    s.potl(n, m, s.r_to_rb(r)) * func(r, l, m);
	else
	  coef[n] += del * qe.weight(i) * pow(r, (int)(s.get_dof()-1)) *
	    s.potl(n, l, s.r_to_rb(r)) * func(r, l, m);
	break;

      case potential:
	if (s.get_dof()==2)
	  coef[n] += del * qe.weight(i) * pow(r, (int)(s.get_dof()-1)) *
	    s.dens(n, m, s.r_to_rb(r)) * func(r, l, m);
	else
	  coef[n] += del * qe.weight(i) * pow(r, (int)(s.get_dof()-1)) *
	    s.dens(n, l, s.r_to_rb(r)) * func(r, l, m);
	break;

      }
    }
  }

  for (int n=0; n<numc; n++)
    coef[n] /= sqrt(s.norm(n, l));

  return coef;

}
