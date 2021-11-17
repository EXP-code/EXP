#include <localmpi.H>
#include <EllipsoidForce.H>
#include <ZBrent.H>

static double find_fct(double u, vector<double>& z)
{
  double ans = -1.0;
  for (int k=0; k<3; k++) ans += z[k]*z[k]/(z[3+k]*z[3+k] + u);
  return ans;
}

int    EllipsoidForce::n          = 1000;
double EllipsoidForce::dfrac      = 2.0e-3;
string EllipsoidForce::cache_file = ".ellip_force.cache";

EllipsoidForce::EllipsoidForce(double a0, double a1, double a2, double mass,
			       double rmin, double rmax,
			       BarType bartype, bool quadr,
			       double param, int num, int numr)
{
  a.push_back(a0);
  a.push_back(a1);
  a.push_back(a2);

  this->mass    = mass; 
  this->N       = num;
  this->numr    = numr;
  this->param   = param;
  this->bartype = bartype;
  this->rmin    = rmin;
  this->rmax    = rmax;
  this->tmade   = false;
  this->quadr   = quadr;

  // *M_PI*a[0]*a[1]*a[2] sucked into the leading factor for the 
  // gravitational potential
  //

  switch (bartype) {
  case powerlaw:
    rho0 = (2.0*param + 3.0)*mass/4.0; 
    break;
  case ferrers:
    rho0 = 2.0*exp(lgamma(2.5+param) - lgamma(1.5) - lgamma(1.0+param))*mass/4.0;
    break;
  case expon:
    rho0 = a[0]*a[0]*mass/(4.0*param*param) / 
      ( 1.0 - (1.0 + a[0]/param)*exp(-a[0]/param) ); 
    break;
  }

  gq = std::make_shared<LegeQuad>(N);
  gz = std::make_shared<LegeQuad>(4*N);

  if (quadr) {

    vector<double> x(3, 0);
    lrmin = log(rmin);
    lrmax = log(rmax);
    ldr = (lrmax - lrmin)/n;

    double msum = 0.0;
    double psum = 0.0;
    double r1=0, r2=0, rlst, rcur=0;

    for (int j=0; j<=n; j++) {
      
      x[0] = exp(lrmin + ldr*j);
      
      if (j) {
	
	r1 = r2;
	r2 = x[0];

	rlst = rcur;
	rcur = RhoBar(r2);
	
	msum += 4.0*M_PI/2.0 * ldr * (rlst*pow(r1, 3) + rcur*pow(r2, 3));
	psum += 4.0*M_PI/2.0 * ldr * (rlst*pow(r1, 2) + rcur*pow(r2, 2));
	
      } else {
	
	r1 = r2;
	rcur = RhoBar(x[0]);
      
      }
    
      rr.push_back(x[0]);
      mm.push_back(msum);
      pp.push_back(psum);
      uu.push_back(U22(x[0]));
    }
  }

}

//
// The Destructor (obviously)
//
EllipsoidForce::~EllipsoidForce()
{
  // Nothing
}


//
// Solve for u or lambda
//
double EllipsoidForce::solve(vector<double> x, double m2)
{
  ZBrent< vector<double> > zbrent;

  double r2 = 0.0, m1 = sqrt(m2);
  for (int k=0; k<3; k++) {
    x[k] /= m1;
    r2 += x[k]*x[k];
  }

  // Assume a_0>a_1>a_2 and therefore the largest possible value of
  // u is given by u_{max} \equiv sum x_k^2/m^2 - a_2^2 

  double umax = r2 - a[2]*a[2];

  if (umax<0.0) {
    cerr << "Out of bounds!!! umax=" << umax
	 << " so no solution is possible [this should not happen!]" << endl;
    return 0.0;
  }
  
  x.insert(x.end(), a.begin(), a.end());
  double ans;
  ZBrent< vector<double> > zb;
  ZBrentReturn ret = zb.find(find_fct, x, 0.0, umax, 1.0e-10, ans);
  
  switch (ret) {
  case Bracket:
    cout << "Root is not bracketed" << endl;
    break;

  case Iteration:
    cout << "Number of iterations exceeded" << endl;
    break;

  case Good:
    break;
  }

  return ans;
}

double EllipsoidForce::getDens(vector<double> x)
{
  double m2=0.0, rho=0.0;
  for (int k=0; k<3; k++) m2 += x[k]*x[k]/(a[k]*a[k]);

  if (m2>1.0) return 0.0;

  switch (bartype) {
  case powerlaw:
    rho = rho0*pow(m2, param);
    break;
  case ferrers:
    rho = rho0*pow(1.0 - m2, param);
    break;
  case expon:
    rho =  rho0*exp(-a[0]*sqrt(m2)/param)/sqrt(m2);
    break;
  }

  return rho/(M_PI*a[0]*a[1]*a[2]);
}

void EllipsoidForce::MassInertia(double& M, vector<double>& I)
{
  M = 0.0;
  I = vector<double>(3, 0.0);
  vector<double> z(3);
  double fac;

  for (int i=0; i<N; i++) {
    z[0] = a[0]*gq->knot(i);

    for (int j=0; j<N; j++) {
      z[1] = a[1]*gq->knot(j);

      for (int k=0; k<N; k++) {
	z[2] = a[2]*gq->knot(k);

	fac = gq->weight(i)*gq->weight(j)*gq->weight(k) * getDens(z); 
	M += fac;
	I[0] += fac * (z[1]*z[1] + z[2]*z[2]);
	I[1] += fac * (z[0]*z[0] + z[2]*z[2]);
	I[2] += fac * (z[0]*z[0] + z[1]*z[1]);
      }
    }
  }

  M *= 8.0*a[0]*a[1]*a[2];
  for (int k=0; k<3; k++) I[k] *= 8.0*a[0]*a[1]*a[2];
}


/*
  Internal potential has the form (Chandra pg. 52, eqn. 93)

  V = \pi G a_1 a_2 a_3 \int^\infty_0\frac{du}{\Delta}
  \left[\Psi(1) - \Psi(m^2(u))\right]

  where

  m^2(u) = \sum_{i=1}^3\frac{x_i^2}{a_i^2 + u}

  and

  \Psi(m^2) = \int_1^{m^2} dm^2\,\rho(m^2)

  The external potential has the form (Cha ndra pg. 51, eqn. 89)

  V = \pi G a_1 a_2 a_3 \int^\infty_\lambda\frac{du}{\Delta}
  \left[\Psi(1) - \Psi(m^2(u))\right]
*/
double EllipsoidForce::getPotl(vector<double> x)
{
  double mshell=0.0, ans=0.0, u, d, t, denom, m2;
  double tmin, tmax=0.5*M_PI;

  //
  // Are we inside or outside?
  //
  double ellip = 0.0;
  for (int k=0; k<3; k++) ellip += x[k]*x[k]/(a[k]*a[k]);
  ellip -= 1.0;

  if (ellip<0.0) {		// Inside
    tmin = 0.0;
  } else {			// Outside
    tmin = atan(solve(x, 1.0));
  }

  for (int i=0; i<N; i++) {
    t = tmin + (tmax - tmin)*gq->knot(i);
    u = tan(t);
    d = cos(t);
    d = 1.0/(d*d);
    
    m2 = 0.0;
    denom = 1.0;
    for (int k=0; k<3; k++) {
      m2 += x[k]*x[k]/(a[k]*a[k] + u);
      denom *= a[k]*a[k]+u;
    }

    switch (bartype) {
    case powerlaw:
      if (param==-1)
	mshell = rho0/(param+1.0)*(1.0 - pow(m2, param+1.0));
      else
	mshell = rho0/(param+1.0)*(1.0 - pow(m2, param+1.0));
      break;
    case ferrers:
      mshell = -rho0/(param+1.0)*pow(1.0-m2, param+1.0);
      break;
    case expon:
      mshell = -2.0*rho0*param/a[0]*(exp(-a[0]/param) - exp(-a[0]*sqrt(m2)/param));
      break;
    }

    ans += d*gq->weight(i) * mshell/sqrt(denom);
  }

  return ans*(tmax - tmin);
}

double EllipsoidForce::getSurfDens(double r)
{
  if (!quadr) return 0.0;

  if (r>=a[0]) return 0.0;

  double ans = 0.0, fac = sqrt(1.0 - r*r/(a[0]*a[0]));
  vector<double> z(3);

  z[0] = r;
  z[1] = 0.0;
  for (int k=0; k<4*N; k++) {
    z[2] = a[2]*fac*gz->knot(k);

    ans += a[2]*fac*gz->weight(k) * getDens(z);
  }
  
  return ans * 2.0;
}


const int nphi = 200, ntheta = 100;

double EllipsoidForce::U22(double r) 
{
  if (!quadr) return 0.0;

  const double numfac = 0.25*sqrt(15.0/(2.0*M_PI));
  vector<double> z(3);
  const double dphi = M_PI/nphi;

  double cosx, sinx, phi, ans=0.0;

  if (gt==0) gt = std::make_shared<LegeQuad>(ntheta);

  for (int i=0; i<nphi; i++) {
    phi = dphi*i;
    for (int j=0; j<ntheta; j++) {
      cosx = gt->knot(j);
      sinx = sqrt(1.0 - cosx*cosx);

      z[0] = r*sinx*cos(phi);
      z[1] = r*sinx*sin(phi);
      z[2] = r*cosx;

      ans += gt->weight(j)*dphi * getPotl(z) * sinx*sinx*cos(2.0*phi);
    }
  }

  return ans*numfac*4.0;
}

double EllipsoidForce::RhoBar(double r) 
{
  if (!quadr) return 0.0;

  vector<double> z(3);
  const double dphi = 0.5*M_PI/nphi;

  double cosx, sinx, phi, ans=0.0;

  if (gt==0) gt = std::make_shared<LegeQuad>(ntheta);

  for (int i=0; i<nphi; i++) {
    phi = dphi*i;
    for (int j=0; j<ntheta; j++) {
      cosx = gt->knot(j);
      sinx = sqrt(1.0 - cosx*cosx);

      z[0] = r*sinx*cos(phi);
      z[1] = r*sinx*sin(phi);
      z[2] = r*cosx;

      ans += gt->weight(j)*dphi * getDens(z);
    }
  }

  return ans*8.0/(4.0*M_PI);
}


double EllipsoidForce::getPotl(double r)
{
  if (!quadr) return 0.0;

  if (r<rr.front()) return pp.front();
  if (r>rr.back() ) return -mm.back()/r;

  double lr = log(r);
  int indx = (int)floor((lr - lrmin)/ldr);
  double a = (lrmin + ldr*(indx+1) - lr)/ldr;
  double b = (lr - lrmin - ldr*indx)/ldr;

  return a*pp[indx] + b*pp[indx+1];
}


bool EllipsoidForce::quadpot(double r, double& f, double& fr, double& M)
{
  if (!quadr) return 0.0;

  if (r<rr.front()) M = mm.front();
  if (r>rr.back() ) M = mm.back();

  if (r<rr.front()) {
    f = fr = 0.0;
    return false;
  }

  if (r>=rr.back()) {
    f = uu.back()*pow(rr.back()/r, 3.0);
    fr = -3.0*f/r;
    return false;
  }

  int numt = rr.size();
  double lr = log(r);
  int indx = max<int>(0, min<int>((int)floor((lr - lrmin)/ldr), numt-2));
  int indx2 = max<int>(1, indx);
  double dlr = (lr - (lrmin + ldr*indx2))/ldr;

  double a = (lrmin + ldr*(indx+1) - lr)/ldr;
  double b = (lr - lrmin - ldr*indx)/ldr;

  f = a*uu[indx] + b*uu[indx+1];
  fr = (uu[indx2+1]*(dlr+0.5) - 2.0*uu[indx]*dlr + uu[indx2-1]*(dlr-0.5))/(r*ldr);
  M = a*mm[indx] + b*mm[indx+1];

  return true;
}


void EllipsoidForce::MakeTable(int n1, int n2, int n3)
{
  if (tmade) return;

  ntab.push_back(n1);
  ntab.push_back(n2);
  ntab.push_back(n3);

  int ntot = n1*n2*n3;

  table = vector<double>(4*ntot, 0.0);
  
  rtmin = rmin;
  rtmax = rmax;

  if (rtmin>0.0) {
    tlog = true;
    rtmin = log(rtmin);
    rtmax = log(rtmax);
  }

  for (int i=0; i<3; i++) dX.push_back((rtmax - rtmin)/ntab[i]);

  if (!read_cache()) {

    unsigned cnt = 0;
    vector<double> x(3), x1(3);
    double h;
    int indx;

    for (int i=0; i<n1; i++)  {
      
      for (int j=0; j<n2; j++) {

	for (int k=0; k<n3; k++) {

	  if (cnt++ % numprocs) continue;

	  x[0] = rtmin + dX[0]*i;
	  if (tlog) x[0] = exp(x[0]);
      
	  x[1] = rtmin + dX[1]*j;
	  if (tlog) x[1] = exp(x[1]);
  
	  x[2] = rtmin + dX[2]*k;
	  if (tlog) x[2] = exp(x[2]);
	
	  indx = (i*n2 + j)*n3 + k;

	  table[4*indx] = getPotl(x);

	  //
	  // 5-point derivative
	  //
	  for (int l=0; l<3; l++) {
	    
	    x1 = x;		// Set the coordinate and step
	    h = dX[l] * dfrac;
	    if (tlog) h *= 2.0*x1[l];
	    
	    x1[l] += h;		// First positive step point
	    table[4*indx+l+1] +=  8.0*getPotl(x1);
	    
	    x1[l] += h;		// Second positive step point
	    table[4*indx+l+1] += -getPotl(x1);
	    
	    x1[l] += -h*3;	// First negative step point
	    table[4*indx+l+1] += -8.0*getPotl(x1);

	    x1[l] += -h;	// Second negative step point
	    table[4*indx+l+1] +=  getPotl(x1);

				// Difference denominator
	    table[4*indx+l+1] /=  12.0*h;
	  }
	}
      }
    }

    if (numprocs>1) {
      MPI_Allreduce(MPI_IN_PLACE, &table[0], 4*ntot, MPI_DOUBLE, MPI_SUM,
		    MPI_COMM_WORLD);
    }

    write_cache();

  }
      
  tmade = true;
}

bool EllipsoidForce::nindx(vector<double>& x, vector<int>& n)
{
  for (int i=0; i<3; i++) {
    if (tlog) {
      if (log(x[i])>=rtmax) return false;
      n[i] = static_cast<int>(floor((log(x[i]) - rtmin)/dX[i]));
    } else {
      if (x[i]>=rtmax) return false;
      n[i] = static_cast<int>(floor((x[i] - rtmin)/dX[i]));
    }

				// Sanity bounds
    n[i] = min<int>(n[i], ntab[i]-2);
    n[i] = max<int>(n[i], 0);
  }

  return true;
}

void EllipsoidForce::TableEval(vector<double> x, vector<double>& force)
{
  //
  // Use an array rather than computing the terms on the fly . . .
  //
  static const unsigned char bits[8][3] = {
    {0, 0, 0},
    {0, 0, 1},
    {0, 1, 0},
    {0, 1, 1},
    {1, 0, 0},
    {1, 0, 1},
    {1, 1, 0},
    {1, 1, 1}
  };

  if (!tmade) {
    // User error . . .
    throw "You can't call this routine until MakeTable is called";
  }

  double f;
  int indx;
  vector<double> x1(x), a(3);
  vector<int> n(3), n1(3);

  //
  //  Put lower limit on the grid interpolation
  //
  for (int i=0; i<3; i++) {
    if (tlog)
      x1[i] = max<double>(fabs(x1[i]), exp(rtmin));
    else
      x1[i] = max<double>(fabs(x1[i]), rtmin);
  }

  if ( nindx(x1, n) ) {
    for (int i=0; i<3; i++)
      a[i] = ( rtmin + dX[i]*(n[i]+1) - (tlog ? log(x1[i]) : x1[i]) )/dX[i];
    
				// Zero out for accumulation in loop below
    for (int j=0; j<4; j++) force[j] = 0.0;

				// Compute the interpolation . . . 
				// Each index l is one of the 8 verticies
				// of the enclosing grid cube
    for (int l=0; l<8; l++) {
				// Compute interpolation coefficients 
				// and array index
      n1 = n;
      f = 1.0;
      for (int i=0; i<3; i++) {
	f *= bits[l][i] ? 1.0 - a[i] : a[i];
	n1[i] += bits[l][i] ? 1 : 0;
      }
      indx = (n1[0]*ntab[1] + n1[1])*ntab[2] + n1[2];
      
				// Accumulation of each term for
				// the potential and each force component
      for (int j=0; j<4; j++) force[j] += f*table[4*indx+j];
    }
				// Sign to take care of octant symmetry
				// in force
    for (int j=0; j<3; j++) force[j+1] *= -x[j]/(fabs(x[j])+1.0e-14);
    
  } else {
				// Spherical force and potential
				//
    double r, r2 = 0.0;
    for (int i=0; i<3; i++) r2 += x[i]*x[i];
    r = sqrt(r2);
    force[0] = -mass/r;		// This one is the potential
    for (int i=0; i<3; i++) force[1+i] = - x[i]*mass/(r*r2);
  }
}


void EllipsoidForce::write_cache()
{
  if (myid) return;

  const char name[] = "EllipsoidForce::write_cache(): ";

  ofstream out(cache_file.c_str());
  if (!out) {
    cerr << name << "error opening cache file, "
	 << "file=" << __FILE__ << " line=" << __LINE__ << endl;
    return;
  }

  int j, ntot=1;

  out.write((const char *)&mass,        sizeof(double));
  out.write((const char *)&rtmin,       sizeof(double));
  out.write((const char *)&rtmax,       sizeof(double));
  out.write((const char *)&param,       sizeof(double));
  out.write((const char *)&(j=bartype), sizeof(int));

  for (int i=0; i<3; i++)
    out.write((const char *)&a[i],      sizeof(double));

  for (int i=0; i<3; i++) {
    out.write((const char *)&ntab[i],   sizeof(int));
    ntot *= ntab[i];
  }

  out.write((const char *)&table[0],    4*ntot*sizeof(double));
}

/*
  Print error messages for unusable cache and return
*/

#define ERR1(N, V0, V1)					    \
  {							    \
    if (myid==0)                                            \
    cerr << "EllipsoidForce::read_cache(): " << N	    \
	 << "(" << V0 << ") != " << N << "1(" << V1	    \
	 << "), recomputing and caching the table" << endl; \
    return false;					    \
  }

#define ERR2(N, I, V0, V1)				    \
  {							    \
    if (myid==0)                                            \
    cerr << "EllipsoidForce::read_cache(): " << N	    \
	 << "[" << I << "](" << V0 << ") != "		    \
	 << N << "1[" << I << "](" << V1 << ")"		    \
	 << "), recomputing and caching the table" << endl; \
    return false;					    \
  }

bool EllipsoidForce::read_cache()
{
  const char name[] = "EllipsoidForce::read_cache(): " ;

  ifstream in(cache_file.c_str());
  if (!in) {
    if (myid==0)
      cerr << name << "error opening cache file, myid=" <<  myid 
	   << " file=" << __FILE__ << " line=" << __LINE__ << endl;
    return false;
  }

  int j, ntot=1;

  double mass1, rtmin1, rtmax1, param1;
  vector<double> a1(3);
  vector<int>    ntab1(3);
  
  in.read((char *)&mass1,        sizeof(double));
  in.read((char *)&rtmin1,       sizeof(double));
  in.read((char *)&rtmax1,       sizeof(double));
  in.read((char *)&param1,       sizeof(double));
  in.read((char *)&j,            sizeof(int));

  for (int i=0; i<3; i++)
    in.read((char *)&a1[i],      sizeof(double));

  for (int i=0; i<3; i++) {
    in.read((char *)&ntab1[i],   sizeof(int));
    ntot *= ntab1[i];
  }

  if (mass != mass1)     ERR1("mass", mass, mass1);

  if (rtmin != rtmin1)   ERR1("rtmin", rtmin, rtmin1);

  if (rtmax != rtmax1)   ERR1("rtmax", rtmax, rtmax1);

  if (param != param1)   ERR1("param", param, param1);

  if ((int)bartype != j) ERR1("bartype", bartype, j);

  for (int i=0; i<3; i++) {
    if (a[i] != a1[i])   ERR2("a", i, a[i], a1[i]);
  }

  for (int i=0; i<3; i++) {
    if (ntab[i] != ntab1[i]) ERR2("ntab", i, ntab[i], ntab1[i]);
  }

  table = vector<double>(4*ntot);

  in.read((char *)&table[0], 4*ntot*sizeof(double));

  return true;
}
