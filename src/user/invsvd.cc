#include <iostream>
#include <cstdlib>
#include <cmath>
#include <complex>

#include <CVector.h>

typedef complex<double> zomplex;

extern "C" int zgesvd_(char *jobu, char *jobvt, int *m, int *n, 
	zomplex *a, int *lda, double *s, zomplex *u, 
	int *ldu, zomplex *vt, int *ldvt, zomplex *work, 
	int *lwork, double *rwork, int *info);



/* 
   Solves an nth-order linear system of equations a x = b using
	SVD decomposition.       a, b, and x are previously allocated as
	a[n][n], b[n] and x[n].  Only x is modified in the routine. 
*/
int inverse_svd(CMatrix& a, CMatrix& x, double threshold)
{

// Make copies of a and b so as not to disturb their contents

  int n = a.getnrows();
  int m = a.getncols();
  if (n != m) {
    cerr << "svd_inverse: Expecting square matrix" << endl;
    exit(-1);
  }
				// Transform to LAPACK


  zomplex *A = new zomplex [n*n];
  double *S = new double [n];
  zomplex *U = new zomplex [m*m];
  zomplex *VT = new zomplex [n*n];
  int ldm = 3*m;
  zomplex *MWORK = new zomplex [ldm];
  int ldr = 3*n*(5*n-4);
  double *RWORK = new double [ldr];
  int INFO;
  
  int i, j;

  for (j=0; j<a.getncols(); j++)
    for (i=0; i<a.getnrows(); i++) {
      A[i+j*m].real() = a[i+1][j+1].real();
      A[i+j*m].imag() = a[i+1][j+1].imag();
    }

				// Perform singular value decomposition

  char Ac[] = "A";
  zgesvd_(Ac, Ac, &n, &m, A, &m, S, U, &m, VT, &n, MWORK, &ldm, 
	  RWORK, &INFO);

				// Transform back

  CMatrix u(1, n, 1, n);
  CMatrix vt(1, n, 1, n);
  CMatrix d(1, n, 1, n);
  d.zero();

  for (j=0; j<a.getncols(); j++)
    for (i=0; i<a.getnrows(); i++) {
      u[i+1][j+1].real() = U[i+j*m].real();
      u[i+1][j+1].imag() = U[i+j*m].imag();
      vt[i+1][j+1].real() = VT[i+j*m].real();
      vt[i+1][j+1].imag() = VT[i+j*m].imag();

      if (i==j)	d[i+1][i+1] = S[i];
    }
  
				// TEST

  /*

  CMatrix test = u * d * vt;

  double tmp, maxx = 0.0;
  for (i=1; i<=a.getnrows(); i++)
    for (j=1; j<=a.getncols(); j++) {
      tmp = fabs(test[i][j] - a[i][j]);
      if (tmp > maxx) maxx = tmp;
    }

  cout << "Max deviation: " << maxx << endl;
  cout << "Singular values: " << endl;
  for (i=0; i<a.getnrows(); i++) cout << i << "  " << S[i] << endl;

  */

				// End TEST


				// Perform back-substitution

  int icnt = 0;

  double max = fabs( d[1][1] * threshold );
  for (i=1; i<=a.getnrows(); i++) 
    if (fabs(d[i][i]) > max) 
      d[i][i] = 1.0/d[i][i];
    else {
      d[i][i] = 0.0;
      icnt++;
    }


  CMatrix v = vt.Adjoint();
  CMatrix ut = u.Adjoint();
  x = v * d * ut;

				// TEST

  /*

  CMatrix test2 = a * x;
  
  int k;
  max = 0.0;
  for (i=1; i<=a.getnrows(); i++)
    for (j=1; j<=a.getncols(); j++) {
      if (i==j) 
	tmp = fabs(test2[i][j] - 1.0);
      else
	tmp = fabs(test2[i][j]);

      if (tmp > max) max = tmp;
    }

  cout << "Max deviation from inverse: " << max << endl;

  */

				// End TEST

  delete [] A;
  delete [] S;
  delete [] U;
  delete [] VT;
  delete [] MWORK;
  delete [] RWORK;

  return icnt;
}
