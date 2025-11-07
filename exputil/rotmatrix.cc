#include <complex>
#include <cstdlib>
#include <cmath>


#define odd(x) ( (int)(x/2)*2 != x )
#define sign(x) ( copysign(1.0, x) )

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

  if (odd(abs(n-m))) fac *= sign(cos12) * sign(sin12);

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
  return sqrt( (2.0*ll+1)/(4.0*M_PI*M_PI) ) * exp(
         (double)mm*log(2.0) + 0.5*(lgamma(1.0+ll-mm) - lgamma(1.0+ll+mm)) +
	 lgamma(0.5*(ll+mm+1)) - lgamma(0.5*(ll-mm)+1.0) ) * cos(0.5*M_PI*(ll+mm));
}


std::complex<double> VeeBeta(int l, int l2, int m, double beta)
{
  return pow(std::complex<double>(0.0,1.0), (double)(m-l2))*rot_matrix(l,l2,m,beta)* Ylm01(l, l2);
}
    

#ifdef TEST

#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

main(int argc, char **argv)
{
  int L=2, M=2, m;
  double beta = 45.0;

  char *prog=argv[0];
  
  while (1) {

    int c = getopt(argc, argv, "l:m:b:");
    if (c==-1) break;

    switch (c) {
    case 'l': 
      L = atoi(optarg);
      break;
    case 'm': 
      M = atoi(optarg);
      break;
    case 'b':
      beta = atof(optarg);
      break;
    }
  }

  beta *= M_PI/180.0;

  cout.precision(5);

  double norm = 0.0, element;
  double plgndr(int l, int m, double x);

  std::complex<double> z;

  for (m=-L; m<=L; m++) {
    z = VeeBeta(L, m, M, beta);
    cout << setw(5)  << m 
	 << setw(14) << (element=rot_matrix(L, M, m, beta))
	 << setw(14) << Ylm01(L, m)*plgndr(L, abs(m), cos(beta))
	 << setw(14) << Re(z)
	 << setw(14) << Im(z)
	 << endl;
    norm += element*element;
  }
  cout << "\nNorm=" << sqrt(norm) << endl;
  
}

#endif

#ifdef TEST2

#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "gaussQ.H"

main(int argc, char **argv)
{
  int L=2, M=2, m=2;
  int nint=40;

  char *prog=argv[0];
  
  while (1) {

    int c = getopt(argc, argv, "L:M:m:n:");
    if (c==-1) break;

    switch (c) {
    case 'L': 
      L = atoi(optarg);
      break;
    case 'M': 
      M = atoi(optarg);
      break;
    case 'm': 
      m = atoi(optarg);
      break;
    case 'n':
      nint = atoi(optarg);
      break;
    }
  }


  LegeQuad lq(nint);


  cout.setf(ios::scientific);
  cout.precision(15);

  double beta, ans = 0.0;
  for (int i=0; i<nint; i++) {
    beta = acos(2.0*(lq.knot(i)-0.5));
    ans += 2.0*lq.weight(i)*rot_matrix(L, m, M, beta)*rot_matrix(L, m, M, beta);
  }
    
  cout << "\nAns=" << ans << endl;
  
}

#endif
