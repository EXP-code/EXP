#include <stdio.h>
#include <math.h>
#include <Vector.h>

#include <kevin_complex.h>


#include <clinalg.h>
#include <octave/CmplxSVD.h>


int linear_solve_svd(CMatrix& a, CVector& b, CVector& x, double threshold)

/* Solves an nth-order linear system of equations a x = b using
	SVD decomposition.       a, b, and x are previously allocated as
	a[n][n], b[n] and x[n].  Only x is modified in the routine. */

{
  int *index;

// Make copies of a and b so as not to disturb their contents

  int n = a.getnrows();
  CMatrix m = a;
  x = b;

				// Transform to OCTAVE

  ComplexMatrix m(a.getnrows, a.getncols);

  int i, j;

  for (i=0; i<a.getnrows; i++)
    for (j=0; j<a.getncols; j++) m(i, j) = a[i+1][j+1];


// Perform singular value decomposition and back-substitution

  double d;
  index = new int[n] - 1;
  if (lu_decomp(m, index, d)==1)
    {
      delete [] (index+1);
      return 1;
    }
  lu_backsub(m, index, x);

  delete [] (index+1);
  return 0;
}

