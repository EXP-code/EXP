// C++
#include <unistd.h>
#include <stdlib.h>
#include <values.h>
#include <math.h>
#include <Vector.h>

#include <iostream>

using namespace std;

#ifdef TEST

int SVD(Matrix A, Matrix &U, Matrix &V, Vector &Z);

int
main(int argc, char **argv)
{
  int i, j, k;
  double tmp, tmp2;
  
  int example=0;
  int diag_only=0;

  if (argc>1) example = atoi(*++argv);
  if (argc>2) diag_only=1;

  int nRow, nCol;
  Matrix a, u, v, b;
  Vector z;

  switch (example) {
  case 0:

    nRow = 3;
    nCol = 3;
    a = Matrix (0,nRow-1,0,nCol-1);
    u = Matrix (0,nRow-1,0,nCol-1);
    v = Matrix (0,nCol-1,0,nCol-1);
    b = Matrix (0,nRow-1,0,nCol-1);
    z = Vector (0,nCol-1);	       

    a[0][0]=0; a[0][1]=0; a[0][2]=0; 
    a[1][0]=0; a[1][1]=0; a[1][2]=0; 
    a[2][0]=0; a[2][1]=0; a[2][2]=0; 
    break;

  case 1:

    nRow = 4;
    nCol = 3;
    a = Matrix (0,nRow-1,0,nCol-1);
    u = Matrix (0,nRow-1,0,nCol-1);
    v = Matrix (0,nCol-1,0,nCol-1);
    b = Matrix (0,nRow-1,0,nCol-1);
    z = Vector (0,nCol-1);	       

    a[0][0]=5; a[0][1]=.000001; a[0][2]=1; 
    a[1][0]=6; a[1][1]=.999999; a[1][2]=1; 
    a[2][0]=7; a[2][1]=2.00001; a[2][2]=1; 
    a[3][0]=8; a[3][1]=2.9999;  a[3][2]=1; 
    break;

  case 2:

    nRow = 2;
    nCol = 2;
    a = Matrix (0,nRow-1,0,nCol-1);
    u = Matrix (0,nRow-1,0,nCol-1);
    v = Matrix (0,nCol-1,0,nCol-1);
    b = Matrix (0,nRow-1,0,nCol-1);
    z = Vector (0,nCol-1);	       


    a[0][0]=1; a[0][1]=3;  
    a[1][0]=-4; a[1][1]=3;

    break;

  case 3:

    nRow = 3;
    nCol = 3;
    a = Matrix (0,nRow-1,0,nCol-1);
    u = Matrix (0,nRow-1,0,nCol-1);
    v = Matrix (0,nCol-1,0,nCol-1);
    b = Matrix (0,nRow-1,0,nCol-1);
    z = Vector (0,nCol-1);	       


    a[0][0]=0.1; a[0][1]=1; a[0][2]=1; 
    a[1][0]=0.100001; a[1][1]=1; a[1][2]=1;
    a[2][0]=1; a[2][1]=1; a[2][2]=2;
    break;

  case 4:
    nRow = 3;
    nCol = 3;
    a = Matrix (0,nRow-1,0,nCol-1);
    u = Matrix (0,nRow-1,0,nCol-1);
    v = Matrix (0,nCol-1,0,nCol-1);
    b = Matrix (0,nRow-1,0,nCol-1);
    z = Vector (0,nCol-1);	       


    a[0][0]=0.1; a[0][1]=1; a[0][2]=1; 
    a[1][0]=0.11; a[1][1]=1; a[1][2]=1;
    a[2][0]=1; a[2][1]=1; a[2][2]=2;
    break;

  case 5:
    nRow = 3;
    nCol = 3;
    a = Matrix (0,nRow-1,0,nCol-1);
    u = Matrix (0,nRow-1,0,nCol-1);
    v = Matrix (0,nCol-1,0,nCol-1);
    b = Matrix (0,nRow-1,0,nCol-1);
    z = Vector (0,nCol-1);	       


    a[0][0]=1; a[0][1]=1; a[0][2]=0; 
    a[1][0]=1; a[1][1]=0; a[1][2]=0;
    a[2][0]=0; a[2][1]=0; a[2][2]=0;
    break;

  case 6:

    nRow = 2;
    nCol = 2;
    a = Matrix (1,nRow,1,nCol);
    u = Matrix (1,nRow,1,nCol);
    v = Matrix (1,nCol,1,nCol);
    b = Matrix (1,nRow,1,nCol);
    z = Vector (1,nCol);	       


    a[1][1]=1; a[1][2]=3;  
    a[2][1]=-4; a[2][2]=3;

    break;

  case 7:
    cin >> nRow;
    cin >> nCol;
    a = Matrix (0,nRow-1,0,nCol-1);
    u = Matrix (0,nRow-1,0,nCol-1);
    v = Matrix (0,nCol-1,0,nCol-1);
    b = Matrix (0,nRow-1,0,nCol-1);
    z = Vector (0,nCol-1);	       

    for (i=0; i<nRow; i++) {
      for (j=0; j<nCol; j++)
	cin >> a[i][j];
    }
    break;

  }
  
  if (!SVD(a, u, v, z)) {
    cerr<< "\nMessage from _main: SVD bombed!\n";
    exit(0);
  }

  if (!diag_only) {

    for (i=0; i<nCol; i++)
      if (z[i]>0.0)
	cout << sqrt(z[i]) << " ";
      else
	cout << z[i] << " ";
    cout << endl << endl;
    for (i=0; i<nRow; i++)
      {
	for (j=0; j<nCol; j++)
	  if (z[j]<1.0e-8)
	    cout << u[i][j] << " ";
	  else
	    cout << u[i][j]/sqrt(z[j]) << " ";
	cout << endl;
      }
    cout << endl;
    for (i=0; i<nCol; i++)
      {
	for (j=0; j<nCol; j++)
	  cout << v[i][j] << " ";
	cout << endl;
      }
  }
  
  double max=0.0;
  double amax=0.0;
  if (!diag_only) cout << "\nTest:\n";
  for (i=0; i<nRow; i++) {
    
    for (j=0; j<nCol; j++) {
      
      b[i][j] = 0.0;
      for (k=0; k<nCol; k++)
	b[i][j] += u[i][k]*v[j][k];
      
      if (!diag_only) cout << b[i][j] << " ";
      tmp = fabs(a[i][j]-b[i][j]);
      max = max < tmp ? tmp : max;
      tmp /= (fabs(a[i][j]) + MINDOUBLE);
      amax = amax < tmp ? tmp : amax;
    }
    if (!diag_only) cout << endl;
  }
  cout << "\nUWV':  Rel error: " << amax << "    Abs error: " << max << "\n\n";

  if (!diag_only) cout.form("\nU:\n");
  
  max=0.0;
  amax=0.0;
  for (i=0; i<nCol; i++) {
    for (j=i; j<nCol; j++) {
      
      tmp = 0.0;
      for (k=0; k<nRow; k++)
	if (z[i]<1.0e-8 || z[j]<1.0e-8)
	  tmp += u[k][i]*u[k][j];
	else
	  tmp += u[k][i]*u[k][j]/sqrt(z[i]*z[j]);
      if (!diag_only) cout << tmp << " ";

      if (i==j) {
	tmp2 = fabs(1.0-tmp);
	max = max < tmp2 ? tmp2 : max;
      }
      else {
	tmp2 = fabs(tmp);
	amax = amax < tmp2 ? tmp2 : amax;
      }	  
    }
    if (!diag_only) cout << endl;
  }
  cout << "\nU:  Diag error: " << max << "    Off-diag error: " << amax << "\n\n";

  if (!diag_only) cout << "\nV:\n";
  
  max=0.0;
  amax=0.0;
  for (i=0; i<nCol; i++) {
    for (j=i; j<nCol; j++) {
      
      tmp = 0.0;
      for (k=0; k<nCol; k++)
	tmp += v[i][k]*v[j][k];
      if (!diag_only) cout << tmp << " ";

      if (i==j) {
	tmp2 = fabs(1.0-tmp);
	max = max < tmp2 ? tmp2 : max;
      }
      else {
	tmp2 = fabs(tmp);
	amax = amax < tmp2 ? tmp2 : amax;
      }	  
    }
    if (!diag_only) cout << cout;
  }
  cout << "\nV:  Diag error: " << max << "    Off-diag error: " << amax << "\n\n";

}

#endif

/* This SVD routine is based on pgs 26-32 of "Compact Numerical Methods
   for Computers" by J. C. Nash (1979), used to compute the pseudoinverse.

   Original implementation by Bryant Marks.

   The routine computes the decomposition:

       A = U D V'

   Inputs:

   A = Matrix(0, m, 0, n)

   Returns:

   U = Matrix(0, m, 0, n)	Prefactor

   V = Matrix(0, n, 0, n)	Inverse (transpose) of postfactor

   Z = Vector(0, n)		Square of the diagonal elements of D


   Returns 1 if sucessful, 0 otherwise.

   Modified by MDW (2/14/93) :

        Translation to C++
        To do orginal Givens rotation
        Support for the Matrix class

   Comments:

   This project was motivated by the fact that I was not able to get
   the Numerical Recipes routines to work in all cases; more folks
   than me had this problem.  This one seems to work in cases where
   the NR routine fails.  Also, I did not implement the "enhancement"
   discussed by Nash because it fails in special cases
   (e.g. successive rows of leading zeros) making it less general.  
   This routine could be *very* easily migrated back to standard C.

*/

				/* MACHINE CONSTANT */
				/* Minimum inverse power of two */
// const double XX_EPS = 1.0/DMAXPOWTWO;
// const double XX_EPS = MINDOUBLE;
const double XX_EPS = 1.11078e-16;

int SVD(Matrix &A, Matrix &U, Matrix &V, Vector &Z)
{
  int i, j, k, RotCount, SweepCount, slimit;
  double eps, tol, vt, p, x0, y0, q, r, c0, s0, d1, d2;
  void bomb_SVD(const char* error);

  if (A.getrlow()!=0 || A.getclow()!=0) {
    bomb_SVD("\tMatrix out of bounds\n\tThis routine assumes input Matrix A(0,n,0,n)");
    return 0;
  }

  int nRow = A.getrhigh() + 1;
  int nCol = A.getchigh() + 1;

  eps = XX_EPS;
  tol = nRow*nRow*eps;

  slimit = nCol/4+1;		// Safety feature to truncate
  if (slimit < 7.0)		// sweep iteraction in case of
    slimit = 7;			// exotic failure
  SweepCount = 0;
  RotCount = 1;

  U = A;

  for (i=0; i<nCol; i++)
    for (j=0; j<nCol; j++)
      {
	V[i][j] = 0.0;
	V[i][i] = 1.0;
      }

  while (RotCount != 0 && SweepCount <= slimit)
    {
      RotCount = nCol*(nCol-1)/2;
      SweepCount++;
      for (j=0; j<nCol-1; j++)
	{
	  for (k=j+1; k<nCol; k++)
	    {
	      p = q = r = 0.0;
	      for (i=0; i<nRow; i++)
		{
		  x0 = U[i][j]; y0 = U[i][k];
		  p += x0*y0; q += x0*x0; r += y0*y0;
		}

	      if (q<r) {	// Exchange columns which are out of order
		c0 = 0.0;
		s0 = 1.0;
		
		for (i=0; i<nRow; i++)
		  {
		    d1 = U[i][j]; 
		    d2 = U[i][k];
		    U[i][j] = d1*c0+d2*s0; 
		    U[i][k] = -d1*s0+d2*c0;
		  }
		
		for (i=0; i<nCol; i++)
		  {
		    d1 = V[i][j]; 
		    d2 = V[i][k];
		    V[i][j] = d1*c0+d2*s0; 
		    V[i][k] = -d1*s0+d2*c0;
		  }
		
	      }
	      else
		{
		  if (fabs(q*r) < tol)
		    RotCount--;
		  else if (p*p/(q*r) <= tol)
		    RotCount--;
		  else
		    {
		      q -= r;
		      vt = sqrt(4*p*p+q*q);
		      c0 = sqrt(fabs(0.5*(vt+q)/vt)); 
		      s0 = p/(vt*c0);

		      for (i=0; i<nRow; i++)
			{
			  d1 = U[i][j]; 
			  d2 = U[i][k];
			  U[i][j] = d1*c0+d2*s0; 
			  U[i][k] = -d1*s0+d2*c0;
			}

		      for (i=0; i<nCol; i++)
			{
			  d1 = V[i][j]; 
			  d2 = V[i][k];
			  V[i][j] = d1*c0+d2*s0; 
			  V[i][k] = -d1*s0+d2*c0;
			}
		    }
		}
	    }
	}
				// Sweep complete
#if DEBUG      
      {
	int imax = nCol*(nCol-1)/2;
	cerr << "Sweep = " << SweepCount << " # of rotations skipped = "
	     << imax - RotCount << " of " << imax  << '\n';
      }
#endif
				// New sweep?
    }

				// Compute singular values

  for (j=0; j<nCol; j++)
    {
      q = 0;
      for (i=0; i<nRow; i++) q += U[i][j]*U[i][j];
/*
      q = sqrt(q);
*/
      Z[j] = q;
/*
      if (q>tol)
	{
	  for (i=0; i<nRow; i++) U[i][j] /= q;
	}
*/
    }

  if (SweepCount > slimit) {
    bomb_SVD("Sweep Limit exceeded");
    return 0;
  }
  else
    return 1;
}


void bomb_SVD(const char* error)
{
  cerr << "SVD Error:\n" << error << '\n';
}

