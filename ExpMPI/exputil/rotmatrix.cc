#include <stdlib.h>
#include <math.h>

#include <kevin_complex.h>


// Ref: Tremaine and Weinberg

double rot_matrix(int l, int m, int n, double beta)
{
  

//  if (l<abs(n) || l<abs(m) || beta<0.0 || beta>M_PI) return 0.0;

  if (l<abs(n) || l<abs(m)) return 0.0;

  double fac, ans = 0.0;
  double cos12 = cos(0.5*beta);
  double sin12 = sin(0.5*beta);

  if (cos12==0.0) cos12 = 1.0e-18;
  if (sin12==0.0) sin12 = 1.0e-18;

  fac = exp( 0.5*(lgamma(1.0+l+n) + lgamma(1.0+l-n) +
		  lgamma(1.0+l+m) + lgamma(1.0+l-m) ) );

  if (abs(n-m)%2) fac *= copysign(1.0, cos12) * copysign(1.0, sin12);

  cos12 = log(fabs(cos12));
  sin12 = log(fabs(sin12));

  for (int t=0; 1; t++, fac*=-1.0) {
    if (t+m-n<0) continue;
    if (l-m<t || l+n<t) break;
    ans += fac * exp( cos12*(2.0*(l-t)+n-m) + sin12*(2.0*t+m-n)
		     -(lgamma(1.0+l-m-t) + lgamma(1.0+l+n-t) +
		       lgamma(1.0+t) + lgamma(1.0+t+m-n) ) );
  }

  return ans;
}



/* angular factors . . . (Y_lm)*/

double Ylm01(int ll, int mm)
{
  mm = abs(mm);
  return sqrt( (2.0*ll+1)/(4.0*M_PI*M_PI) ) * exp(
         (double)mm*log(2.0) + 0.5*(lgamma(1.0+ll-mm) - lgamma(1.0+ll+mm)) +
	 lgamma(0.5*(ll+mm+1)) - lgamma(0.5*(ll-mm)+1.0) );
}


Complex VeeBeta(int l, int l2, int m, double beta)
{
  return pow(Complex(0.0,1.0), (double)(m-l2))*rot_matrix(l,l2,m,beta)*
    Ylm01(l, l2);
}
    

#ifdef TEST

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <GetOpt.h>

main(int argc, char **argv)
{
  int L=2, M=2, m;
  double beta = 45.0;

  char c, *prog=argv[0];
  GetOpt getopt(argc, argv, "l:m:b:");

  while ((c = getopt ()) != EOF)
      switch (c) {
      case 'l': 
	L = atoi(getopt.optarg);
	break;
      case 'm': 
	M = atoi(getopt.optarg);
	break;
      case 'b':
	beta = atof(getopt.optarg);
	break;
      }

  beta *= M_PI/180.0;

  cout.precision(5);

  double norm = 0.0, element;
  for (m=-L; m<=L; m++) {
    cout << setw(5) << m << 
      setw(14) << (element=rot_matrix(L, M, m, beta)) << endl;
    norm += element*element;
  }
  cout << "\nNorm=" << sqrt(norm) << endl;
  
}

#endif
