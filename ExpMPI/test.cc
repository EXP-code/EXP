#include <stdlib.h>
#include <math.h>
#include <stl.h>

#include <iostream>
#include <iomanip>
#include <strstream>
#include <fstream>

#include <Interp.h>

struct Pair {
  double r;
  double m;
} cur;

class comp {
public:
  bool operator() (const Pair& x, const Pair& y) const
    {
      return x.r < y.r;;
    }
};

int main()
{
  vector<Pair> pairs;
  char cbuf[256];

  ifstream in("file.new");
  in.getline(cbuf, 256);

  int num;
  {
    istrstream istr(cbuf);
    istr >> num;
  }

  int c;
  double m, x, y, z, r;
  for (int i=0; i<num; i++) {
    in.getline(cbuf, 256);
    istrstream istr(cbuf);
    istr >> c;
    istr >> m;
    istr >> x;
    istr >> y;
    istr >> z;
 
    cur.r = sqrt(x*x + y*y + z*z);
    cur.m = m;
    pairs.push_back(cur);
 }
    
  sort(pairs.begin(), pairs.end(), comp());
 

  vector<double> rt(pairs.size()); // Create radius vector
  vector<double> mt(pairs.size()); // Create cumulated mass vector
  double cum = 0.0;
  for (int i=0; i<pairs.size(); i++) {
    rt[i] = pairs[i].r;
    mt[i] = cum + 0.5*pairs[i].m;
    cum += pairs[i].m;
  }

  int N = pairs.size()/50;
  if (N>1000) N = 1000;

  vector<double> rr(N); // Create radius vector
  vector<double> mm(N); // Create cumulated mass vector
  vector<double> m2(N); // 2nd derivative vector
  vector<double> rho(N); // Create cumulated mass vector
  vector<double> wrk(N); // Create cumulated mass vector
  vector<double> wrk2(N); // Create cumulated mass vector
  vector<double> pot(N); // Create cumulated mass vector

  double RMIN = rt[rt.size()-1]*0.001;
  double dy, dr = (rt[rt.size()-1] - rt[0] - RMIN)/N;
  for (int i=0; i<N; i++) {
    rr[i] = rr[0] + dr*i + RMIN;
    mm[i] = odd2(rr[i], rt, mt);
  }

  Spline(rr, mm, 0.0, 0.0, m2);

  for (int i=0; i<N; i++) {
    Splint2(rr, mm, m2, rr[i], y, dy);
    rho[i] = dy/(4.0*M_PI*rr[i]*rr[i]);
    wrk[i] = -dy/rr[i];
  }

  Spline(rr, wrk, 0.0, 0.0, wrk2);
  splsum2(rr, wrk, wrk2, pot);
  for (int i=0; i<N; i++) wrk[i] = pot[i];
  for (int i=0; i<N; i++) pot[i] = wrk[N-1] - wrk[i];
  for (int i=1; i<N; i++) pot[i] += -mm[i]/rr[i];

		// Test
  dr = (rt[rt.size()-1] - rt[0])/N;
  for (int i=0; i<N; i++) {
    r = rr[0] + dr*i;
    Splint2(rr, mm, m2, r, y, dy);

    cout << setw(15) << r
	 << setw(15) << y
	 << setw(15) << dy
	 << setw(15) << odd2(r, rr, rho)
	 << setw(15) << odd2(r, rr, pot)
	 << endl;
  }
}
