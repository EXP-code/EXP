#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include <EllipForce.H>

std::string EllipForce::cache_file = ".ellipforce.cache";

EllipForce::EllipForce(double A, double B, double C, double MASS,
		       int NUM, int NUMR)
{
  a = A;
  b = B;
  c = C;
  num = NUM;
  numr = NUMR;
  mass = MASS;

  if (read_cache()) return;

  lq = std::make_shared<LegeQuad>(num);

  double x, y, z, dr = a/(numr-1);
  double ans, xfac, yfac, zfac;

  // Work vectors
  std::vector<double> dm(numr), w1(numr), w2(numr);

  // Radius, Mass, Potential
  r = std::vector<double>(numr);
  m = std::vector<double>(numr);
  p = std::vector<double>(numr);

  double mfac = mass/(4.0*M_PI/3.0*a*b*c);

  for (int v=0; v<numr; v++) {

    r[v] = dr*v;

    ans = 0.0;
    xfac = min<double>(r[v], a);
    for (int i=1; i<=num; i++) {
      x = xfac*lq->knot(i);
      for (int j=1; j<=num; j++) {
	yfac = sqrt(xfac*xfac - x*x);
	y = yfac*lq->knot(j);
	for (int k=1; k<=num; k++) {
	  zfac = sqrt(xfac*xfac - x*x - y*y);
	  z = zfac*lq->knot(k);

	  // Is inside ellipsoid?
	  if (x*x/(a*a) + y*y/(b*b) + z*z/(c*c) < 1.0) 
	    ans += xfac*yfac*zfac*lq->weight(i)*lq->weight(j)*lq->weight(k);
	}
      }
    }
    
    m[v] = 8.0*ans*mfac;	// Tabulate total mass
  }

				// External potential integrand: (dM/dr)/r
  for (int v=0; v<numr; v++) {
    if (r[v] <= 0.0) w1[v] = 0.0;
    else w1[v] = drv2(r[v], r, m)/r[v];
  }

				// Integrate external potential 
  w2[0] = 0.0;			// using trapezoidal rule
  for (int v=1; v<numr; v++) 
    w2[v] = 0.5*(r[v] - r[v-1])*(w1[v] + w1[v-1]) + w2[v-1];

				// Compute the total gravitational potential
  for (int v=0; v<numr; v++) {
    if (r[v] <= 0.0) p[v] = -w2[numr-1];
    else p[v] = -m[v]/r[v] - (w2[numr-1] - w2[v]);
  }

  write_cache();

}
  
void EllipForce::write_cache()
{
  const char efmagic[] = "efcache";
  std::ofstream out(cache_file.c_str());
  if (!out) {
    cerr << "EllipForce: could not write cache file" << endl;
    return;
  }

  out.write((const char *)efmagic, 7*sizeof(char));
  out.write((const char *)&a,    sizeof(double));
  out.write((const char *)&b,    sizeof(double));
  out.write((const char *)&c,    sizeof(double));
  out.write((const char *)&num,  sizeof(int));
  out.write((const char *)&numr, sizeof(int));
  out.write((const char *)&mass, sizeof(double));

  out.write((const char *)&r[0], numr*sizeof(double));
  out.write((const char *)&m[0], numr*sizeof(double));
  out.write((const char *)&p[0], numr*sizeof(double));
}

bool EllipForce::read_cache()
{
  const char efmagic[] = "efcache", efmagicT[] = "efcache";

  ifstream in(cache_file.c_str());
  if (!in) {
    cerr << "EllipForce: could not read cache file" << endl;
    return false;
  }

				// Check magic
  in.read((char *)efmagicT, 7*sizeof(char));
  for (int i=0; i<7; i++)
    if (efmagic[i] != efmagicT[i]) return false;

				// Check parameters
  double A, B, C, MASS;
  int NUM, NUMR;
  in.read((char *)&A,    sizeof(double));
  in.read((char *)&B,    sizeof(double));
  in.read((char *)&C,    sizeof(double));
  in.read((char *)&NUM,  sizeof(int));
  in.read((char *)&NUMR, sizeof(int));
  in.read((char *)&MASS, sizeof(double));

  if (A   != a   || B   != b    || C    != c || 
      NUM != num || NUMR!= numr || MASS != mass) {
    cout << "EllipForce: cache has different parameters" << endl;
    return false;
  }

				// Read arrays
  r = std::vector<double>(numr);
  m = std::vector<double>(numr);
  p = std::vector<double>(numr);

  in.read(( char *)&r[0], numr*sizeof(double));
  in.read(( char *)&m[0], numr*sizeof(double));
  in.read(( char *)&p[0], numr*sizeof(double));

				// Check for good stream
  if (in) {
    cout << "EllipForce: cache read" << endl;
    return true;
  }
  else {
    cout << "EllipForce: bad cache" << endl;
    return false;
  }
}

EllipForce::~EllipForce()
{
  // NADA
}

double EllipForce::getMass(double x)
{
  if (x>=r[numr-1]) return m[numr-1];
  return odd2(x, r, m);
}

double EllipForce::getPot(double x)
{
  if (x>a) return -m[numr-1]/x;
  return odd2(x, r, p);
}

void EllipForce::PrintTable()
{
				// Print the table
  for (int v=0; v<numr; v++)
    cout << setw(15) << r[v]
	 << setw(15) << m[v]
	 << setw(15) << p[v]
	 << endl;
}

void EllipForce::PrintModel(const string& file)
{
  ofstream out(file.c_str());
  if (!out) return;

  out << numr << endl;
				// Compute density
				// 
  std::vector<double> d(numr);
  for (int v=1; v<numr; v++) {
    d[v] = drv2(r[v], r, m);
    d[v] /= 4.0*M_PI*r[v]*r[v];
  }
				// Extrapolate to get inner point
				// 
  d[0] = (d[1]*(r[2] - r[0]) + d[2]*(r[0] - r[1]))/(r[2] - r[1]);

  for (int v=0; v<numr; v++) {
    out << setw(15) << r[v]
	<< setw(15) << d[v]
	<< setw(15) << m[v]
	<< setw(15) << p[v]
	<< endl;
  }
}

void EllipForce::TestTable()
{
				// Print the table
  for (int v=0; v<numr; v++)
    cout << setw(15) << r[v]
	 << setw(15) << m[v] - getMass(r[v])
	 << setw(15) << p[v] - getPot(r[v])
	 << endl;
}
