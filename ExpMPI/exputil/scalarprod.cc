// This may look like C code, but it is really -*- C++ -*-

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


#include <math.h>
#include <Vector.h>
#include <gaussQ.h>
#include <biorth2d.h>
#include <OrthoPoly.h>


Vector scalar_prod(ScalarType type, double rmin, double rmax, int l, int m,
		   AxiSymBiorth& s, double (*func)(double, int, int), 
		   int numc, int numg)
{

  if ( !(s.get_dof()==2 || s.get_dof()==3) ) {
    cerr << "scalar_prod: dof=" << s.get_dof() << ",  must be 2 or 3!\n";
    exit(-1);
  }

  LegeQuad qe(numg);

  Vector coef(1, numc);
  double r, del = rmax - rmin;
  coef.zero();

  int i, n;

  for (i=1; i<=numg; i++) {
    r = rmin + del*qe.knot(i);
    for (n=1; n<=numc; n++) {
      switch (type) {

      case density:
	if (s.get_dof()==2)
	  coef[n] += del * qe.weight(i) * pow(r, (int)(s.get_dof()-1)) *
	    s.potl(n, m, s.r_to_rb(r)) * (*func)(r, l, m);
	else
	  coef[n] += del * qe.weight(i) * pow(r, (int)(s.get_dof()-1)) *
	    s.potl(n, l, s.r_to_rb(r)) * (*func)(r, l, m);
	break;

      case potential:
	if (s.get_dof()==2)
	  coef[n] += del * qe.weight(i) * pow(r, (int)(s.get_dof()-1)) *
	    s.dens(n, m, s.r_to_rb(r)) * (*func)(r, l, m);
	else
	  coef[n] += del * qe.weight(i) * pow(r, (int)(s.get_dof()-1)) *
	    s.dens(n, l, s.r_to_rb(r)) * (*func)(r, l, m);
	break;

      }
    }
  }

  for (n=1; n<=numc; n++)
    coef[n] /= sqrt(s.norm(n-1, l));

  return coef;

}
