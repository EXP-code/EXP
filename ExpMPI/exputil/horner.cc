/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine returns the results of a Horner's scheme for a Poly
 *  Ref: Henrici, Applied and Computational Complex Analysis, p. 433f
 *
 *  Call sequence:
 *  -------------
 *
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *  Value
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 01/07/93
 *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <Vector.h>
#include <kevin_complex.h>

//#include "poly.h"
#include "cpoly.h"

Vector get_horner(double z, Poly &p)
{
  int i, j;
  int order = p.getorder();

  Matrix b(-1,order,0,order);
  for (i=0; i<=order; i++) b[-1][i] = p[order-i];
  
  for (i=0; i<=order; i++) {
    b[i][0] = b[-1][0];
    for (j=1; j<=order-i; j++)
      b[i][j] = b[i-1][j] + z*b[i][j-1];
  }

  Vector ans(0,order);
  for (i=0; i<=order; i++) ans[i] = b[i][order-i];
  
  return ans;
}


CVector Cget_horner(Complex z, CPoly &p)
{
  int i, j;
  int order = p.getorder();
  CVector ans;

  if (order==0) {
    ans = CVector(0,1);
    ans[0] = p[0];
    ans[1] = 0.0;
    return ans;
  }

  CMatrix b(-1,order,0,order);
  for (i=0; i<=order; i++) b[-1][i] = p[order-i];
  
  for (i=0; i<=order; i++) {
    b[i][0] = b[-1][0];
    for (j=1; j<=order-i; j++)
      b[i][j] = b[i-1][j] + z*b[i][j-1];
  }

  ans = CVector(0,order);
  for (i=0; i<=order; i++) ans[i] = b[i][order-i];
  
  return ans;
}

