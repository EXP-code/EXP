#include <EllipForce.H>

EllipForce::EllipForce(double A, double B, double C, double MASS,
		       int NUM, int NUMR)
{
  a = A;
  b = B;
  c = C;
  num = NUM;
  numr = NUMR;
  mass = MASS;

  lq = std::make_shared<LegeQuad>(num);

  double x, y, z, dr = a/numr;
  double ans, xfac, yfac, zfac;

  // Work vectors
  vector<double> dm(numr), w1(numr), w2(numr);

  // Radius, Mass, Potential
  r = vector<double>(numr);
  m = vector<double>(numr);
  p = vector<double>(numr);

  double mfac = mass/(4.0*M_PI/3.0*a*b*c);

  for (int v=0; v<numr; v++) {

    r[v] = dr*v;

    ans = 0.0;
    xfac = min<double>(r[v], a);
    for (int i=0; i<num; i++) {
      x = xfac*lq->knot(i);
      for (int j=0; j<num; j++) {
	yfac = sqrt(xfac*xfac - x*x);
	y = yfac*lq->knot(j);
	for (int k=0; k<num; k++) {
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

}
  
EllipForce::~EllipForce()
{
  // Nothing
}

double EllipForce::getMass(double x)
{
  if (x>a) return m[numr-1];
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

void EllipForce::TestTable()
{
				// Print the table
  for (int v=0; v<numr; v++)
    cout << setw(15) << r[v]
	 << setw(15) << m[v] - getMass(r[v])
	 << setw(15) << p[v] - getPot(r[v])
	 << endl;
}
