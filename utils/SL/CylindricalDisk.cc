#include <CylindricalDisk.H>

#include <sstream>

extern double plgndr(int, int, double);


CylindricalDisk::~CylindricalDisk()
{
  // NADA
}

void CylindricalDisk::Initialize(double Rmin, double Rmax, bool logR,
				 int Nmax, int Lmax,
				 int numR, int numT, int numG,
				 vector<double> param)
{
  if (param.size()) SetParameters(param);

  //
  // First compute spherical model from cylindrical distribution
  //
  LegeQuad lr(numR), lt(numT);
  vector<double> r(numR), d(numR), m(numR,0), p(numR), p2(numR,0);
  vector<double> X(3, 0);
  double trmin = 0.3*Rmin;
  double trmax = Rmax;
    
  if (logR) {
    trmin = log(0.3*Rmin);
    trmax = log(Rmax);
  }

  double dRR = (trmax - trmin)/(numR - 1);

  // Compute density
  for (int i=0; i<numR; i++) {
    r[i] = trmin + dRR*i;
    if (logR) r[i] = exp(r[i]);
    d[i] = 0.5*Density(r[i])/r[i];
  }
  
  // Compute mass and potential integral
  for (int i=1; i<numR; i++) {
    double dr = r[i] - r[i-1];
    m[i]  = m[i-1]  + 2.0*M_PI*dr*(r[i]*r[i]*d[i] + r[i-1]*r[i-1]*d[i-1]);
    p2[i] = p2[i-1] + 2.0*M_PI*dr*(r[i]*d[i] + r[i-1]*d[i-1]);
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
  numt = numT;
  numg = numG;
  
  model = std::make_shared<SphericalModelTable>
    (numR, 
     &r[0]-1, &d[0]-1, &m[0]-1, &p[0]-1,
     0, 0, 0, "Shells from cylinder");
  
  ortho = std::make_shared<SLSphere>
    (lmax, nmax, numr, rmin, rmax, 1, 1, model);

  // Compute coefficients

  dRR = log(rmax) - log(rmin);

  coefs = vector< Eigen::VectorXd >(lmax+1);
  for (int l=0; l<=lmax; l++) {
    coefs[l].resize(nmax);
    coefs[l].setZero();
  }
  
  double theta, rho;
  double x, y, z, R;
  Eigen::VectorXd vec(nmax);

  for (int i=0; i<numr; i++) {

    R = rmin*exp(dRR*lr.knot(i));
    
    for (int k=1; k<=numt; k++) {
	
      theta = acos(2.0*(lt.knot(k)-0.5));
	
      x = X[0] = R * sin(theta);
      y = X[1] = 0.0;
      z = X[2] = R * cos(theta);
	
      R = sqrt(x*x + y*y + z*z);

      //
      // 4Pi needed for correct sol'n to Poisson's eqn
      //
      rho = 4.0*M_PI * Density(X);
      
      for (int l=0; l<=lmax; l++) {

	//
	// Get ortho fcts
	//
	ortho->potl(nmax, l, R, vec);

	coefs[l] += 2.0*lt.weight(k) * dRR*R*R*R*lr.weight(i) * 
	  2.0*M_PI * (0.5*l+0.25)/M_PI * plgndr(l, 0, 2.0*(lt.knot(k)-0.5)) *
	  rho * vec;
      }
    }
  }

  initialized = true;
}

void CylindricalDisk::make_grids()
{
  if (!initialized) 
    throw "CylindricalDisk::make_grids(): you must call Initialize(...) before attemping any field evaluations!";

  linear = (rmin<=0.0) ? true : false;
  if (linear) {
    rgmin = rmin;
    rgmax = rmax;
  } else {
    rgmin = log(rmin);
    rgmax = log(rmax);
  }
  
  //
  // Use 3% mesh spacing for derivative stencil
  //
  const double frac = 0.03;
  dR = (rgmax - rgmin)/(numg - 1);

  frc = vector< vector< vector<double> > >(numg);
  for (int i=0; i<numg; i++) {
    frc[i] = vector< vector<double> >(numg);
    for (int j=0; j<numg; j++)
      frc[i][j] = vector<double>(4);
  }
  
  //
  // Compute derivatives with 5-pt formula
  //
  double r, z, hr, hz;
  for (int i=0; i<numg; i++) {
    for (int j=0; j<numg; j++) {
      r = rgmin + dR*i; 
      z = rgmin + dR*j;
      if (!linear) {
	r = exp(r);
	z = exp(z);
      }
      
      hr = r*frac;
      hz = z*frac;

      frc[i][j][0] = density_eval(r, 0.0, z);

      frc[i][j][1] = potential_eval(r, 0.0, z);

      frc[i][j][2] = -(
		          -potential_eval(r+2.0*hr, 0.0, z)
		      +8.0*potential_eval(r+hr    , 0.0, z)
		      -8.0*potential_eval(r-hr    , 0.0, z)
		          +potential_eval(r-2.0*hr, 0.0, z)
		      )/(12.0*hr);

      frc[i][j][3] = -(
		          -potential_eval(r, 0.0, z+2.0*hz)
		      +8.0*potential_eval(r, 0.0, z+hz    )
		      -8.0*potential_eval(r, 0.0, z-hz    )
			  +potential_eval(r, 0.0, z-2.0*hz)
		      )/(12.0*hz);
    }
  }
  
  indx = vector<int>   (4);
  aint = vector<double>(4);
  grid = true;
}

double CylindricalDisk::density_eval(double x, double y, double z)
{
  double r     = sqrt(x*x + y*y + z*z);
  double costh = z/(r+1.0e-18);
  double ans   = 0.0;

  for (int l=0; l<=lmax; l++) {
    ans += -ortho->get_dens(r, l, coefs[l]) * 
      plgndr(l, 0, costh) / (4.0*M_PI);
  }
  
  return ans;
}


double CylindricalDisk::potential_eval(double x, double y, double z)
{
  double r     = sqrt(x*x + y*y + z*z);
  double costh = z/(r+1.0e-18);
  double ans   = 0.0;

  for (int l=0; l<=lmax; l++) {
    ans += ortho->get_potl(r, l, coefs[l]) * plgndr(l, 0, costh);
  }
    
  return ans;
}

void CylindricalDisk::force_eval(double x, double y, double z, 
				 double &fr, double &fz)
{
  if (!grid) make_grids();

  double rr = sqrt(x*x + y*y);
  double zz = fabs(z);

  if (get_interp(rr, zz)) {
    fr = 
      aint[0] * aint[2] * frc[indx[0]][indx[2]][2] +
      aint[1] * aint[2] * frc[indx[1]][indx[2]][2] +
      aint[0] * aint[3] * frc[indx[0]][indx[3]][2] +
      aint[1] * aint[3] * frc[indx[1]][indx[3]][2] ;

    fz = 
      aint[0] * aint[2] * frc[indx[0]][indx[2]][3] +
      aint[1] * aint[2] * frc[indx[1]][indx[2]][3] +
      aint[0] * aint[3] * frc[indx[0]][indx[3]][3] +
      aint[1] * aint[3] * frc[indx[1]][indx[3]][3] ;

    if (z<0.0) fz *= -1.0;

  } else {
    double r3 = sqrt(rr*rr + zz*zz);
    double ff = -model->get_mass(r3)/(r3*r3*r3);
    fr = -ff*rr;
    fz = -ff*z;
  }
}


void CylindricalDisk::pot_dens_eval(double x, double y, double z, 
				    double &potl, double &dens)
{
  if (!grid) make_grids();

  double rr = sqrt(x*x + y*y);
  double zz = fabs(z);

  if (get_interp(rr, zz)) {

    dens = 
      aint[0] * aint[2] * frc[indx[0]][indx[2]][0] +
      aint[1] * aint[2] * frc[indx[1]][indx[2]][0] +
      aint[0] * aint[3] * frc[indx[0]][indx[3]][0] +
      aint[1] * aint[3] * frc[indx[1]][indx[3]][0] ;

    potl = 
      aint[0] * aint[2] * frc[indx[0]][indx[2]][1] +
      aint[1] * aint[2] * frc[indx[1]][indx[2]][1] +
      aint[0] * aint[3] * frc[indx[0]][indx[3]][1] +
      aint[1] * aint[3] * frc[indx[1]][indx[3]][1] ;

  } else {
    double r3 = sqrt(rr*rr + zz*zz);
    dens = 0.0;
    potl = -model->get_mass(r3)/r3;
  }
}


bool CylindricalDisk::get_interp(double x, double z)
{
  if (!linear) {
    if (x==0.0) x = rgmin;
    else        x = log(x);

    if (z==0.0) z = rgmin;
    else        z = log(z);
  }

  if (x>=rgmax || z>=rgmax) return false;

  if (x<rgmin) {
    indx[0] = 0;
    indx[1] = 1;
    aint[0] = 1.0;
    aint[1] = 0.0;
  } else {
    indx[0] = static_cast<int>(floor((x - rgmin)/dR));
    indx[1] = indx[0] + 1;
    aint[0] = (rgmin + dR*indx[1] - x)/dR;
    aint[1] = (x - rgmin - dR*indx[0])/dR;
  }

  if (z<rgmin) {
    indx[2] = 0;
    indx[3] = 1;
    aint[2] = 1.0;
    aint[3] = 0.0;
  } else {
    indx[2] = static_cast<int>(floor((z - rgmin)/dR));
    indx[3] = indx[2] + 1;
    aint[2] = (rgmin + dR*indx[3] - z)/dR;
    aint[3] = (z - rgmin - dR*indx[2])/dR;
  }

  return true;
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
