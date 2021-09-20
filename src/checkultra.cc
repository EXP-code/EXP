#include <math.h>
#include <iostream>
#include <iomanip>

void usage(char *dummy){};

double ultra(int, double, double);
void get_ultra(int, double, double, Vector& p);

int main()
{
  double x, tmp, del;
  int nmax, l;

  del = 0.0001;
  x = 0.35;
  l = 4;
  nmax = 12;

  Eigen::VectorXd u(nmax);
  Eigen::VectorXd du(nmax);

  get_ultra(nmax, (double)l,   x, u);
  get_ultra(nmax, (double)(l+1),   x, du);
  for (i=nmax; i>0; i--) du[i] = du[i-1];
  du[0] = 0.0;

  for (int i=0; i<nmax; i++) {
    tmp = (ultra(i, l, x+del) - ultra(i, l, x-del))/(2.0*del);
    cout << setw(4) << i << setw(15) << tmp << setw(15) << 
      2.0*(l+1)*du[i-1] << "\n";
  }

}
