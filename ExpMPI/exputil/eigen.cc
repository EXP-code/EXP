using namespace std;

#include <unistd.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <numerical.h>
#include <Vector.h>

static void bomb_ghql(const char* msg)
{
  cerr << msg << endl;
  exit(-1);
}

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double dpythag(double a, double b)
{
  double absa, absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}

static void tred2(Matrix& a, int n, Vector& d, Vector& e)
{
  int l, k, j, i;
  double scale, hh, h, g, f;
  
  for (i=n;i>=2;i--) {
    l=i-1;
    h = scale = 0.0;
    if (l > 1) {
      for (k=1;k<=l;k++)
	scale += fabs(a[i][k]);
      if (scale == 0.0)
	e[i] = a[i][l];
      else {
	for (k=1;k<=l;k++) {
	  a[i][k] /= scale;
	  h += a[i][k]*a[i][k];
	}
	f = a[i][l];
	g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
	e[i] = scale*g;
	h -= f*g;
	a[i][l] = f-g;
	f = 0.0;
	for (j=1;j<=l;j++) {
	  a[j][i] = a[i][j]/h;
	  g = 0.0;
	  for (k=1;k<=j;k++)
	    g += a[j][k]*a[i][k];
	  for (k=j+1;k<=l;k++)
	    g += a[k][j]*a[i][k];
	  e[j] = g/h;
	  f += e[j]*a[i][j];
	}
	hh = f/(h+h);
	for (j=1;j<=l;j++) {
	  f = a[i][j];
	  e[j] = g = e[j]-hh*f;
	  for (k=1;k<=j;k++)
	    a[j][k] -= (f*e[k]+g*a[i][k]);
	}
      }
    } else
      e[i] = a[i][l];
    d[i] = h;
  }
  d[1] = 0.0;
  e[1] = 0.0;
  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[i]=a[i][i]; */
  for (i=1;i<=n;i++) {
    l = i-1;
    if (d[i]) {
      for (j=1;j<=l;j++) {
	g = 0.0;
	for (k=1;k<=l;k++)
	  g += a[i][k]*a[k][j];
	for (k=1;k<=l;k++)
	  a[k][j] -= g*a[k][i];
      }
    }
    d[i] = a[i][i];
    a[i][i] = 1.0;
    for (j=1;j<=l;j++) a[j][i] = a[i][j] = 0.0;
  }
}


static void tqli(Vector& d, Vector& e, int n, Matrix& z)
{
  double pythag(double a, double b);
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;
  
  for (i=2;i<=n;i++) e[i-1]=e[i];
  e[n]=0.0;
  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
	dd=fabs(d[m])+fabs(d[m+1]);
	if ((double)(fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
	if (iter++ == 30) bomb_ghql("tqli: too many iterations");
	g=(d[l+1]-d[l])/(2.0*e[l]);
	r=pythag(g, 1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r, g));
	s=c=1.0;
	p=0.0;
	for (i=m-1;i>=l;i--) {
	  f=s*e[i];
	  b=c*e[i];
	  e[i+1]=(r=pythag(f, g));
	  if (r == 0.0) {
	    d[i+1] -= p;
	    e[m]=0.0;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  d[i+1]=g+(p=s*r);
	  g=c*r-b;
	  for (k=1;k<=n;k++) {
	    f=z[k][i+1];
	    z[k][i+1]=s*z[k][i]+c*f;
	    z[k][i]=c*z[k][i]-s*f;
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l]=g;
	e[m]=0.0;
      }
    } while (m != l);
  }
}


static void eigsrt(Vector& d, Matrix& v, int n)
{
  int k, j, i;
  double p;
  
  for (i=1;i<n;i++) {
    p=d[k=i];
    for (j=i+1;j<=n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=1;j<=n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}


Vector Matrix::Symmetric_Eigenvalues(Matrix &EigenVec) throw (char const*)
{
  double **v, *d, **a;
  int n, niters;
  
  if (rlow != 1 || clow != 1)
    {
      puts("Eigenvalues expects a unit-offset matrix");
      exit(0);
    }
  
  if (rhigh != chigh)
    {
      puts("I can't take the eigenvalues of a non-square matrix");
      exit(0);
    }
  
  n = rhigh;
  
  a = nr_matrix(1, n, 1, n);
  v = nr_matrix(1, n, 1, n);
  d = nr_vector(1, n);
  
  
  
  
  // copy the input matrix into a double pointer array
  
  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
      a[i][j] = (*this)[i][j];
  
  // call numerical recipes eigenroutines
  
  jacobi(a, n, d, v, &niters);
  eigsrt(d, v, n);
  
  
  // copy the eigenvector array into a matrix object
  
  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
      EigenVec[i][j] = v[i][j];
  
  // copy the eigenvalue return into a vector object
  
  
  Vector EigenVals(1, n, d);
  
  
  
  // tidy up memory
  
  free_nr_matrix(a, 1, n, 1, n);
  free_nr_matrix(v, 1, n, 1, n);
  free_nr_vector(d, 1, n);
  
  
  // return
  
  return EigenVals;
}


Vector Symmetric_Eigenvalues(Matrix &m, Matrix &ev) throw (char const*)
{
  Vector tmp;
  
  tmp = m.Symmetric_Eigenvalues(ev);
  
  return tmp;
}


Vector Matrix::Symmetric_Eigenvalues_GHQL(Matrix &EigenVec) throw (char const*)
{
  int n;
  
  if (rlow != 1 || clow != 1)
    {
      puts("Eigenvalues expects a unit-offset matrix");
      exit(0);
    }
  
  if (rhigh != chigh)
    {
      puts("I can't take the eigenvalues of a non-square matrix");
      exit(0);
    }
  
  n = rhigh;
  
  Matrix a(1, n, 1, n);
  Vector d(1, n);
  Vector e(1, n);
  
  a = (*this);
  
  tred2(a, n, d, e);
  tqli(d, e, n, a);
  eigsrt(d, a, n);
  
  // copy the eigenvectors return into a return matrix
  EigenVec = a;
  
  // return
  
  return d;
}


Vector Symmetric_Eigenvalues_GHQL(Matrix &m, Matrix &ev) throw (char const*)
{
  Vector tmp;
  
  tmp = m.Symmetric_Eigenvalues(ev);
  
  return tmp;
}
