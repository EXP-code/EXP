#include <CylindricalDisk.H>

#include <sstream>

extern double plgndr(int, int, double);


CylindricalDisk::~CylindricalDisk()
{
  delete model;
  delete ortho;
}

CylindricalDisk::CylindricalDisk(double Rmin, double Rmax, int Lmax,
				 int Nmax, int numR, int numt)
{
  //
  // First compute spherical model from cylindrical distribution
  //
  LegeQuad lq(numt);
  double dR = (Rmax - 0.3*Rmin)/(numR - 1), cost, sint;
  vector<double> r(numR), d(numR), m(numR,0), p(numR), p2(numR,0);
  vector<double> X(3, 0);

  // Compute density
  for (int i=0; i<numR; i++) {
    r[i] = 0.3*Rmin + dR*i;
    d[i] = 0.0;
    for (int t=1; t<=numt; t++) {
      cost = lq.knot(t);
      sint = sqrt(1.0 - cost*cost);
      X[0] = sint*r[i];
      X[2] = cost*r[i];
      d[i] += lq.weight(t) * Density(X);
    }
  }
  
  // Compute mass and potential integral
  for (int i=1; i<numR; i++) {
    m[i]  = m[i-1]  + 2.0*M_PI*dR*(r[i]*r[i]*d[i] + r[i-1]*r[i-1]*d[i-1]);
    p2[i] = p2[i-1] + 2.0*M_PI*dR*(r[i]*d[i] + r[i-1]*d[i-1]);
  }

  // Compute potential
  for (int i=0; i<numR; i++) {
    if (r[i]<=0.0)
      p[i] = p2[i] - p2[numR-1];
    else
      p[i] = -m[i]/r[i] + p2[i] - p2[numR-1];
  }

  // Debug
  const string tfile("disk_sphere.model");
  ofstream out(tfile.c_str());
  if (out) {
  } else {
    cerr << "Could not open <" << tfile << ">" << endl;
  }

  out << "#" << endl;
  out.precision(10);
  for (int i=0; i<numR; i++) {
    out << setw(18) << r[i]
	<< setw(18) << d[i]
	<< setw(18) << m[i]
	<< setw(18) << p[i]
	<< endl;
  }
  // End Debug

  rmin = Rmin;
  rmax = Rmax;
  lmax = Lmax;
  nmax = Nmax;
  numr = numR;

  model = new SphericalModelTable(numR, 
				  &r[0]-1, &d[0]-1, &m[0]-1, &p[0]-1,
				  0, 0, 0, "Shells from cylinder");
  
  ortho = new SLSphere(lmax, nmax, numr, rmin, rmax, 1, 1, model);


  // Compute coefficients

  dR = (log(rmax) - log(rmin))/numr;

  coefs = vector< Vector >(lmax+1);
  for (int l=0; l<=lmax; l++) coefs[l] = Vector(1, nmax);
  
  double theta, phi, rho, facp;
  double x, y, z, R;
  Vector vec(1, nmax);

  for (int i=0; i<numr; i++) {

    R = rmin*exp(dR*i);
    
    for (int k=1; k<=numt; k++) {
	
      theta = acos(2.0*(lq.knot(k)-0.5));
	
      x = X[0] = R * sin(theta);
      y = X[1] = 0.0;
      z = X[2] = R * cos(theta);
	
      R = sqrt(x*x + y*y + z*z);

      // 4Pi needed for correct sol'n to Poisson's eqn
      rho = 4.0 * M_PI * Density(X);
      
      for (int l=0; l<=lmax; l++) {

	// Get ortho fcts
	ortho->potl(nmax, l, R, vec);

	coefs[l] += 2.0*lq.weight(k) *
	  (0.5*l+0.25)/M_PI *
	  plgndr(l, 0, 2.0*(lq.knot(k)-0.5)) * rho * vec;
      }
    }
  }
}

double CylindricalDisk::density_eval(double x, double y, double z)
{
  double r = sqrt(x*x + y*y + z*z);
  double costh = z/(r+1.0e-18);
  double phi = atan2(y, x);

  double ans = 0.0;
  for (int l=0; l<=lmax; l++) {
    ans += -ortho->get_dens(r, l, coefs[l]) * 
      plgndr(l, 0, costh) / (4.0*M_PI);
  }
  
  return ans;
}


double CylindricalDisk::potential_eval(double x, double y, double z)
{
  double r = sqrt(x*x + y*y + z*z);
  double costh = z/(r+1.0e-18);
  double phi = atan2(y, x);

  double ans = 0.0;
  for (int l=0; l<=lmax; l++) {
    ans += ortho->get_potl(r, l, coefs[l]) * plgndr(l, 0, costh);
  }
    
  return ans;
}

void CylindricalDisk::dump_coefficients(string file)
{
  ofstream out(file.c_str());
  if (!out) {
    cerr << "CylindricalDisk::dump_coefficients: can't open <"
	 << file << "> . . . " << endl;
    return;
  }
  out.setf(ios::scientific);
  out.precision(8);

  double fac;

  //
  // Header
  //
  out << "# Axisymmetric (m=0) coefficients only" << endl;
  out << "#" << setw(4) << "l";
  for (int i=1; i<=nmax; i++) {
    ostringstream ostr;
    ostr << "n=" << i;
    out << setw(18) << ostr.str();
  }
  out << endl;

  //
  // Coefficients
  //
  for (int l=0; l<=lmax; l++) {
    out << setw(5) << l;
    
    fac = sqrt((0.5*l+0.25)/M_PI);
    for (int i=1; i<=nmax; i++) out << setw(18) << coefs[l][i]/fac;
      
    out << endl;
  }
}
