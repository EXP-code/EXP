/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  This routine . . .
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
 *  MDW 11/20/91
 *
 ***************************************************************************/

#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <Vector.h>
#include <kevin_complex.h>

/*
*	Computes integral using Simpson's rule with Romberg's
*	correction with a supplied integrand routine
*
*	Synopsis:
*		answer = rombe2(a,b,f,n);
*               Complex answer, f();
*		double	a,b,f();
*		int	n;
*
* 	Use:
* 		a	initial pt.
* 		b	final pt.
* 		f	function
* 		ans	contains answer on return
* 		t	holds first column of Romberg tableux
* 
*/
#define         DMACH 1.0e-8
#define         NMAX  20

static Complex tc[NMAX];
static double  tr[NMAX];

Complex Crombe(double a, double b, CVector& f)
{
  Complex               ends, x2, x4;
  double         	xmesh,d;
  int			nmax,num,i,j,k,index;

  int n1 = f.getlow();
  int n2 = f.gethigh();
  int n = (int)( log((double)(n2-n1))/log(2.0) + DMACH );
  nmax = (int)(pow(2.0,(double)n) + DMACH);
  if (n2-n1 != nmax) {
    cerr << "Crombe: wrong # of points\n";
    exit (-1);
  }
  if (n > NMAX) {
    cerr << "Crombe: dimension too large\n";
    exit (-1);
  }

  xmesh = (b-a)/nmax;

  d = (b-a)/2.0;
  ends = f[n1] + f[n2];
  num = nmax/2;
  x2 = f[n1+num];
  tc[0] = d*(ends+4.0*x2)/3.0;
  for(i=2; i<=n; i++)
    {
      x4 = 0.0;
      d = d/2.0;
      index = num/2;
      for(j=1; j<=nmax/num; j++)
	{
	  x4 = x4 + f[n1 +index];
	  index = index + num;
	}
      tc[i-1] = d*(ends+2.0*x2+4.0*x4)/3.0;
      
      x2 = x2+x4;
      num = num/2;
      
    }
  
  num=3;
  for(j=1; j <= n-1; j++)
    {
      for(k=1; k <= n-j; ++k)
	tc[k-1] = tc[k]+(tc[k]-tc[k-1])/num;
      num = 4*num+3;
    }

  return tc[0];
}

double Vrombe(double a, double b, Vector& f)
{
  double		ends, x2, x4;
  double         	xmesh,d;
  int			nmax,num,i,j,k,index;

  int n1 = f.getlow();
  int n2 = f.gethigh();
  int n = (int)( log((double)(n2-n1))/log(2.0) + DMACH );
  nmax = (int)(pow(2.0,(double)n) + DMACH);
  if (n2-n1 != nmax) {
    cerr << "Crombe: wrong # of points\n";
    exit (-1);
  }
  if (n > NMAX) {
    cerr << "Vrombe: dimension too large\n";
    exit (-1);
  }

  xmesh = (b-a)/nmax;

  d = (b-a)/2.0;
  ends = f[n1] + f[n2];
  num = nmax/2;
  x2 = f[n1+num];
  tr[0] = d*(ends+4.0*x2)/3.0;
  for(i=2; i<=n; i++)
    {
      x4 = 0.0;
      d = d/2.0;
      index = num/2;
      for(j=1; j<=nmax/num; j++)
	{
	  x4 = x4 + f[n1 +index];
	  index = index + num;
	}
      tr[i-1] = d*(ends+2.0*x2+4.0*x4)/3.0;
      
      x2 = x2+x4;
      num = num/2;
      
    }
  
  num=3;
  for(j=1; j <= n-1; j++)
    {
      for(k=1; k <= n-j; ++k)
	tr[k-1] = tr[k]+(tr[k]-tr[k-1])/num;
      num = 4*num+3;
    }

  return tr[0];
}

#undef          DMACH
