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

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>

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

#include <Eigen/Eigen>

static std::complex<double> tc[NMAX];
static double tr[NMAX];

std::complex<double> Crombe(double a, double b, Eigen::VectorXcd& f)
{
  std::complex<double> ends, x2, x4;

  int sz = f.size();
  int n  = (int)( log((double)(sz-1))/log(2.0) + DMACH );
  int nmax = (int)(pow(2.0,(double)n) + DMACH);
  if (sz != nmax) {
    std::cerr << "Crombe: wrong # of points\n";
    exit (-1);
  }
  if (n > NMAX) {
    std::cerr << "Crombe: dimension too large" << std::endl;
    exit (-1);
  }

  double d = (b-a)/2.0;
  ends = f[0] + f[sz-1];
  int num = nmax/2;
  x2 = f[num];
  tc[0] = d*(ends+4.0*x2)/3.0;
  for (int i=2; i<=n; i++)
    {
      x4 = 0.0;
      d = d/2.0;
      int index = num/2;
      for (int j=1; j<=nmax/num; j++)
	{
	  x4 = x4 + f[index];
	  index = index + num;
	}
      tc[i-1] = d*(ends+2.0*x2+4.0*x4)/3.0;
      
      x2 = x2+x4;
      num = num/2;
    }
  
  num=3;
  for (int j=1; j <= n-1; j++)
    {
      for(int k=1; k <= n-j; ++k)
	tc[k-1] = tc[k]+(tc[k]-tc[k-1])/static_cast<double>(num);
      num = 4*num+3;
    }

  return tc[0];
}

double Vrombe(double a, double b, Eigen::VectorXf& f)
{
  double ends, x2, x4;

  int sz = f.size();
  int n  = (int)( log((double)(sz-1))/log(2.0) + DMACH );
  int nmax = (int)(pow(2.0,(double)n) + DMACH);
  if (sz-1 != nmax) {
    std::cerr << "Crombe: wrong # of points" << std::endl;
    exit (-1);
  }
  if (n > NMAX) {
    std::cerr << "Vrombe: dimension too large" << std::endl;
    exit (-1);
  }

  double d = (b-a)/2.0;
  ends = f[0] + f[sz-1];
  int num = nmax/2;
  x2 = f[num];
  tr[0] = d*(ends+4.0*x2)/3.0;
  for (int i=2; i<=n; i++)
    {
      x4 = 0.0;
      d = d/2.0;
      int index = num/2;
      for(int j=1; j<=nmax/num; j++)
	{
	  x4 = x4 + f[index];
	  index = index + num;
	}
      tr[i-1] = d*(ends+2.0*x2+4.0*x4)/3.0;
      
      x2 = x2+x4;
      num = num/2;
      
    }
  
  num=3;
  for(int j=1; j <= n-1; j++)
    {
      for(int k=1; k <= n-j; ++k)
	tr[k-1] = tr[k]+(tr[k]-tr[k-1])/num;
      num = 4*num+3;
    }

  return tr[0];
}

#undef          DMACH
