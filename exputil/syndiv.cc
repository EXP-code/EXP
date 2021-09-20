/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Synthetic division routines
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
#include "cpoly.H"

void syn_poly_div(Poly &p1, Poly &p2, Poly &Quotient, Poly &Remainder)
{
  int k, j;
  int n1 = p1.getorder();
  int n2 = p2.getorder();
  
  Poly r = p1;
  Poly q = Poly(n1);

  for (k=n1-n2; k>=0; k--) {
    q[k] = r[n2+k]/p2[n2];
    for (j=n2+k-1; j>=k; j--)
      r[j] -= q[k]*p2[j-k];
  }
  
  for (j=n2;j<=n1;j++) r[j]=0.0;

  Quotient = q;
  Remainder = r;
}


void Csyn_poly_div(CPoly &p1, CPoly &p2, CPoly &Quotient, CPoly &Remainder)
{
  int k, j;
  int n1 = p1.getorder();
  int n2 = p2.getorder();
  
  CPoly r(p1);
  CPoly q(n1);

  for (k=n1-n2; k>=0; k--) {
    q[k] = r[n2+k]/p2[n2];
    for (j=n2+k-1; j>=k; j--)
      r[j] -= q[k]*p2[j-k];
  }
  
  for (j=n2;j<=n1;j++) r[j]=0.0;

  Quotient = q;
  Remainder = r;
}


