#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <limits>

#include <global.H>
#include <ResPotOrb.H>
#include <ZBrent.H>
#include <localmpi.H>

#include <pthread.h>  
static pthread_mutex_t iolock = PTHREAD_MUTEX_INITIALIZER;

#undef DEBUG_VERBOSE		// More copious reporting on errors
#undef DEBUG_DEBUG		// Only for a few-step diagnostic run

double ResPotOrb::ALPHA = 0.25;
double ResPotOrb::DELTA_E = 0.001;
double ResPotOrb::DELTA_K = 0.001;
double ResPotOrb::DELTA_B = 0.001;
double ResPotOrb::DELE = 0.001;
double ResPotOrb::DELK = 0.001;
double ResPotOrb::TOLITR = 1.0e-8;
int ResPotOrb::NUME = 200;
int ResPotOrb::NUMX = 200;
int ResPotOrb::RECS = 100;
int ResPotOrb::ITMAX = 50;
std::complex<double> ResPotOrb::I(0.0, 1.0);

const char* ResPotOrb::ReturnDesc[] = {
  "OK", 
  "CoordRad", "CoordVel", "CoordCirc", "CoordTP", "CoordRange", 
  "CoordBadPos", "CoordBadVel", 
  "UpdateBadL", "UpdateIterate", "UpdateBadVal"};


ResPotOrb::ResPotOrb(AxiSymModPtr mod, double mass,
		     int l, int m, int l1, int l2)
{
  halo_model = mod;
  
  grid_computed = false;
  second_order = true;
  dbgFILE = "update.dbg";
  
  L = l;
  M = m;
  L1 = l1;
  L2 = l2;
  
  MASS = mass;

  Rmin = halo_model->get_min_radius();
  Rmax = halo_model->get_max_radius();
  Emax = halo_model->get_pot(Rmax)*(1.0+DELTA_E);
  Emin = halo_model->get_pot(Rmin)*(1.0-DELTA_E);
  Kmin = DELTA_K;
  Kmax = 1.0-DELTA_K;
  Kupd = 0.0;
  
  // SphericalOrbit::RMAXF=1.0;
  orb = std::make_shared<SphericalOrbit>(halo_model);
  orb->set_numerical_params(RECS);
  
  compute_grid();
}

ResPotOrb::~ResPotOrb()
{
  // Nothing
}

static AxiSymModel *hmod;

// Circular orbit
//
static double get_rc(double r, double& J)
{
  return J*J - hmod->get_mass(r)*r;
}

// Turning points
//
static double adj_r(double r, vector<double>& p)
{
  return p[1] - 0.5*p[0]*p[0]/(r*r) - hmod->get_pot(r);
}


void ResPotOrb::compute_grid() 
{
  if (grid_computed) return;
  
  double E, K;
  
  // Ang mom of minimum radius circular orbit
  minJ = Rmin*sqrt(max<double>(halo_model->get_dpot(Rmin)*Rmin, 0.0));
  
  // Ang mom of maximum radius circular orbit
  maxJ = Rmax*sqrt(max<double>(halo_model->get_dpot(Rmax)*Rmax, 0.0));
  
  // For derivatives
  delE = (Emax - Emin)*0.5*DELE;
  delK = (Kmax - Kmin)*0.5*DELK;
  
  double J, rcirc, Ecirc;
  
  dX = 1.0/(NUMX-1);
  
  hmod = halo_model.get();
  
  int iret;
  ZBrent<double> circ;
  
  for (int s=0; s<NUMX; s++) {
    
    // Unit interval to J
    J = Jx(double(s)/(NUMX-1), minJ, maxJ);
    
    // Compute Ecirc
    iret = circ.find(get_rc, J, Rmin, Rmax, 1.0e-10, rcirc);
    if (iret) {
      cerr << "Error locating circular orbit!!\n";
      exit(-1);
    }
    Ecirc = 0.5*J*J/(rcirc*rcirc) + halo_model->get_pot(rcirc);
    if (Ecirc<Emin) Ecirc = Emin;
    
    EminX.push_back(Ecirc);
    
    vector<double> IJ1, EJ1;
    ovector orbvec;
    
    // Grid from Ecirc to Emax
    double dE = (Emax - Ecirc)/(NUME - 1);
    
    for (int i=0; i<NUME; i++) {
      
      E = Ecirc + dE*i;
      
      orb->new_orbit(E, 0.5);
      K = J/orb->Jmax();
      K = max<double>(K, Kmin);
      K = min<double>(K, Kmax);
      
      RW rw;
      
      orb->new_orbit(E, K);
      struct ANGLE_GRID * grid = orb->get_angle_grid();
      LegeQuad Gkn(grid->num);
      
      for (int j=0; j<grid->num; j++) {
	rw.ff.push_back(
			Gkn.weight(j+1)*grid->dw1dt(0, j) *
			cos(grid->w1(0, j)*L1 + grid->f(0, j)*L2));
	rw.w1.push_back(grid->w1(0, j));
	rw.f.push_back(grid->f(0, j));
	rw.r.push_back(grid->r(0, j));
      }

      rw.num = grid->num;
      rw.O1 = orb->get_freq(0);
      rw.O2 = orb->get_freq(1);
      rw.I1 = orb->get_action(0);
      rw.E = E;
      rw.K = K;
      rw.Jm = orb->Jmax();
      
      EJ1.push_back(rw.E);
      IJ1.push_back(rw.I1);
      
      
      // Derivatives
      // -----------
      
      double Ep=E+delE, Em=E-delE, Kp=K+delK, Km=K-delK;
      
      // Energy bounds
      if (E+delE > Emax) {
	Ep = Emax;
	Em = Emax - delE*2.0;
      }
      
      if (K-delK < Kmin) {
	Kp = Kmin + delK*2.0 ;
	Km = Kmin;
      }
      
      // Kappa bounds
      if (K+delK > Kmax) {
	Kp = Kmax;
	Km = Kmax - delK*2.0;
      }
      
      if (K-delK < Kmin) {
	Kp = Kmin + delK*2.0 ;
	Km = Kmin;
      }
      
      // Energy deriv
      orb->new_orbit(Ep, K);
      grid = orb->get_angle_grid();
      for (int j=0; j<grid->num; j++) {
	rw.ff_Ep.push_back(
			   Gkn.weight(j+1)*grid->dw1dt(0, j)*
			   cos(grid->w1(0, j)*L1 + grid->f(0, j)*L2));
	rw.r_Ep.push_back(grid->r(0, j));
      }
      rw.dJm = orb->Jmax();
      
      orb->new_orbit(Em, K);
      grid = orb->get_angle_grid();
      for (int j=0; j<grid->num; j++) {
	rw.ff_Em.push_back(
			   Gkn.weight(j+1)*grid->dw1dt(0, j)*
			   cos(grid->w1(0, j)*L1 + grid->f(0, j)*L2));
	rw.r_Em.push_back(grid->r(0, j));
      }
      
      // Kappa deriv
      orb->new_orbit(E, Kp);
      grid = orb->get_angle_grid();
      for (int j=0; j<grid->num; j++) {
	rw.ff_Kp.push_back(
			   Gkn.weight(j+1)*grid->dw1dt(0, j)*
			   cos(grid->w1(0, j)*L1 + grid->f(0, j)*L2));
	rw.r_Kp.push_back(grid->r(0, j));
      }

      orb->new_orbit(E, Km);
      grid = orb->get_angle_grid();
      for (int j=0; j<grid->num; j++) {
	rw.ff_Km.push_back(
			   Gkn.weight(j+1)*grid->dw1dt(0, j)*
			   cos(grid->w1(0, j)*L1 + grid->f(0, j)*L2));
	rw.r_Km.push_back(grid->r(0, j));
      }
      
      // Finally, store the data for this phase space point
      orbvec.push_back(rw);
      check_rw(rw);
    }
    
    // Normalize energy to [0,1]
    for (size_t j=0; j<EJ1.size(); j++) 
      EJ1[j] = (EJ1[j] - Ecirc)/(Emax - Ecirc);
    
    EX.push_back(EJ1);
    I1X.push_back(IJ1);
    
    orbmat.push_back(orbvec);
  }
  
  numx = EX.size();
  ngrid = orb->get_angle_grid()->num;
  
  grid_computed = true;
}





ResPotOrb::ReturnCode ResPotOrb::coord(double* ps1, double* vel,
				 double& E, double& K, double& I1, double& J,
				 double& O1, double& O2,
				 double& W1, double& W2, double& W3, 
				       double& F, double& BETA, double& PSI)
{
  double pos[3];
  for (int i=0; i<3; i++) pos[i] = ps1[i];
  if (fabs(ps1[2])<1.0e-8) pos[2] = 1.0e-8;
  
  // Compute polar coords
  // --------------------
  
  double r2=0.0, v2=0.0, rv=0.0, r;
  for (int i=0; i<3; i++) {
    r2 += pos[i]*pos[i];
    v2 += vel[i]*vel[i];
    rv += pos[i]*vel[i];
  }
  
  r = sqrt(r2);
  
  if (r>halo_model->get_max_radius()) return CoordRad;
  
  double theta = acos(pos[2]/r);
  double phi = atan2(pos[1], pos[0]);
  
  // Compute E, J, beta
  // ------------------
  
  E = 0.5*v2 + halo_model->get_pot(r);
  
  double angmom[3];
  
  angmom[0] = pos[1]*vel[2] - pos[2]*vel[1];
  angmom[1] = pos[2]*vel[0] - pos[0]*vel[2];
  angmom[2] = pos[0]*vel[1] - pos[1]*vel[0];
  
  J = 0.0;
  for (int i=0; i<3; i++) J += angmom[i]*angmom[i];
  J = sqrt(J);
  
  BETA = 0.0;
  if (J>0.0) BETA = acos(angmom[2]/J);
  
  
  if (fabs(rv)>1.0e-10*sqrt(r2*v2)) rv /= sqrt(r2*v2);
  else {			// Apocenter
    if (J*J/r2 < halo_model->get_mass(r))
      rv = -1.0;
    else			// Pericenter
      rv = 1.0;
  }
  
  // Linear interpolation coefficients
  // ---------------------------------
  
  J = max<double>(J, minJ);
  J = min<double>(J, maxJ);
  
  double X = xJ(J, minJ, maxJ);
  
  int indxX = (int)( X/dX );
  indxX = min<int>(indxX, numx-2);
  indxX = max<int>(indxX, 0);
  
  double cX[2];
  cX[0] = (dX*(indxX+1) - X)/dX;
  cX[1] = 1.0 - cX[0];
  
  double Y, Ecirc = cX[0]*EminX[indxX] + cX[1]*EminX[indxX+1];
  
  if (E<Ecirc) Y = 0.0;
  else if (E>Emax) Y = 1.0;
  else Y = (E - Ecirc)/(Emax - Ecirc);
  
  int indxE;
  double cE[2];
  for (int i1=0; i1<2; i1++) {
    indxE = Vlocate(Y, EX[indxX+i1]);
    indxE = min<int>(indxE, EX[indxX+i1].size()-2);
    indxE = max<int>(indxE, 0);
    
    cE[0] = (EX[indxX+i1][indxE+1] - Y)/
      (EX[indxX+i1][indxE+1] - EX[indxX+i1][indxE]);
    
    // Bounds: keep it on the grid else unphsical!
    cE[0] = max<double>(cE[0], 0.0);
    cE[0] = min<double>(cE[0], 1.0);
    
    cE[1] = 1.0 - cE[0];
    
  }
  
  
  // Compute angles and transform
  // ----------------------------
  
  double fac;
  int num;
  
  vector<double> tF(ngrid, 0.0), tw(ngrid, 0.0), tf(ngrid, 0.0);
  vector<double> tF_Ep(ngrid, 0.0), tF_Em(ngrid, 0.0);
  vector<double> tF_Kp(ngrid, 0.0), tF_Km(ngrid, 0.0);

  vector<double> tR(ngrid, 0.0);
  vector<double> tR_Ep(ngrid, 0.0), tR_Em(ngrid, 0.0);
  vector<double> tR_Kp(ngrid, 0.0), tR_Km(ngrid, 0.0);
  
  I1 = K = O1 = O2 = 0.0;
  
  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<2; i2++) {
      
      RW *rw = &(orbmat[indxX+i1][indxE+i2]);
      num = rw->num;
      fac = cX[i1]*cE[i2];
      
      if (ngrid != num) {
	cerr << "ResPotOrb::coord[1]: Oops! ngrid=" 
	     << ngrid << "  num=" << num 
	     << "  indxX=" << indxX+i1 << "/" << numx
	     << "  indxE=" << indxE+i2 << "/" << EX[indxX+i1].size()
	     << endl;
      }
      
      I1 += fac * rw->I1;
      K  += fac * rw->K;
      O1 += fac * rw->O1;
      O2 += fac * rw->O2;
      
      for (int k=0; k<ngrid; k++) {
	tw[k] += fac * rw->w1[k];
	tf[k] += fac * rw->f[k];
	tF[k] += fac * rw->ff[k];
	tR[k] += fac * rw->r[k];
      }
      
    }
  }
  
  
  // Angle computation
  // -----------------
  
  double psi, w3 = atan2(angmom[1], angmom[0]) + 0.5*M_PI;
  
  if (fabs(BETA)<1.0e-10) psi = phi - w3;
  else {
    double tmp = cos(theta)/sin(BETA);
    
    if (fabs(tmp)>1.0) {
      if (tmp>1.0) psi =  0.5*M_PI;
      else         psi = -0.5*M_PI;
    } 
    else psi = asin(tmp);
  }
  
  
  // W1 and f(W1)
  // ------------

  double w1, f;			// Enforce limits
  if (r < tR[0]) {
    w1 = tw[0];
    f  = tf[0];
  } else if (r > tR[ngrid-1]) {
    w1 = tw[ngrid-1];
    f  = tf[ngrid-1];
  } else {
    w1 = odd2(r, tR, tw);
    f  = odd2(r, tR, tf);
  }
  
  if (rv < 0.0) {
    w1 = 2.0*M_PI - w1;
    f *= -1.0;
  }
  
  // Map Psi into [-Pi, Pi]
  // ----------------------
  double tmp = atan2(sin(phi - w3), cos(phi - w3));
  if (tmp>0.5*M_PI || tmp<-0.5*M_PI) {
    psi = M_PI - psi;
  }
  
  double w2 = psi + f;
  
  K = max<double>(K, Kmin);
  K = min<double>(K, Kmax);
  
  W1  = w1;
  W2  = w2;
  W3  = w3;
  F   = f;
  PSI = psi;
  
  return OK;
}


//
// Coordindate mapping to unit interval
//

double ResPotOrb::xJ(double J, double Jmin, double Jmax)
{
  if (J>=Jmax) return 1.0;
  if (J<=Jmin) return 0.0;
  return pow( (J-Jmin)/(Jmax-Jmin), ALPHA );
}

double ResPotOrb::Jx(double x, double Jmin, double Jmax)
{
  if (x>=1.0) return Jmax;
  if (x<=0.0) return Jmin;
  return Jmin + (Jmax - Jmin)*pow(x, 1.0/ALPHA);
}


double ResPotOrb::dxJ(double J, double Jmin, double Jmax)
{
  if (J>=Jmax) return 1.0/(Jmax-Jmin);
  if (J<=0.0) return std::numeric_limits<float>::max();
  return pow( (J-Jmin)/(Jmax-Jmin), ALPHA-1.0 )/(Jmax - Jmin);
}


void ResPotOrb::getInterp(double I1, double I2, int& indxX, int& indxE,
		       double cX[2], double cE[2], bool& noboundary)
{
  noboundary = true;
  
  // Linear interpolation coefficients
  // ---------------------------------
  
  double X = xJ(I2, minJ, maxJ);
  
  indxX = (int)( X/dX );
  indxX = min<int>(indxX, numx-2);
  indxX = max<int>(indxX, 0);
  
  cX[0] = (dX*(indxX+1) - X)/dX;
  cX[1] = 1.0 - cX[0];
  
  if (std::isnan(cX[0]) || std::isnan(cX[1])) {
    cerr << "getInterp: cX is NaN, X=" << X 
	 << "  I1=" << I1
	 << "  I2=" << I2 
	 << "  minJ=" << minJ
	 << "  maxJ=" << maxJ
	 << endl;
  }
  
  double I1min = cX[0]*I1X[indxX].front()  + cX[1]*I1X[indxX+1].front();
  double I1max = cX[0]*I1X[indxX].back()   + cX[1]*I1X[indxX+1].back();
  
  // Assign Y index
  //
  if (I1<=I1min) {		// Below grid
    noboundary = false;
    indxE = 0;
    cE[0] = 1.0;
    cE[1] = 0.0;
  }
  else if (I1>=I1max) {		// Above grid
    noboundary = false;
    indxE = NUME-2;
    cE[0] = 0.0;
    cE[1] = 1.0;
  }
  else {			// Do binary search to get Y index
    int mid, lo=-1, hi=NUME;
    while ( hi-lo > 1 ) {
      mid=(hi+lo) >> 1;		// Divide by two
      if ( I1 > cX[0]*I1X[indxX][mid] + cX[1]*I1X[indxX+1][mid])
	lo = mid;		// Discard lower interval
      else
	hi = mid;		// Discard upper interval
    }
    
    if (lo==-1)
      indxE = 0;
    else
      indxE = lo;
    
    double Ilo = cX[0]*I1X[indxX][indxE]   + cX[1]*I1X[indxX+1][indxE];
    double Ihi = cX[0]*I1X[indxX][indxE+1] + cX[1]*I1X[indxX+1][indxE+1];
    
    cE[0] = (Ihi - I1)/(Ihi - Ilo);
    cE[1] = 1.0 - cE[0];
    
    // Test: should always be on grid
    //
    if (cE[0]<0.0 && cE[0]>1.0) {
      cout << "ResPotOrb[id=" << myid 
	   << "]: WEIGHT out of bounds, indxE=" << indxE
	   << " cE[0]=" << cE[0]
	   << " cE[1]=" << cE[1]
	   << " Ilo, Ihi, I1=" << Ilo << ", " << Ihi << ", " << I1 << endl;
    }
    if (indxE<0 && indxE>NUME-1) {
      cout << "ResPotOrb[id=" << myid 
	   << "]: INDEX out of bounds, indxE=" << indxE
	   << " cE[0]=" << cE[0]
	   << " cE[1]=" << cE[1]
	   << " Ilo, Ihi, I1=" << Ilo << ", " << Ihi << ", " << I1 << endl;
    }
  }
  
  return;
}

bool ResPotOrb::getValues(double I1, double I2,
			  double& O1, double& O2)
{
  O1 = 0.0;
  O2 = 0.0;
  
  int indxX;
  int indxE;
  double cX[2];
  double cE[2];
  
  bool noboundary;
  getInterp(I1, I2, indxX, indxE, cX, cE, noboundary);
  
  // Interpolate frequencies
  // -----------------------
  
  double fac;
  int num;
  
  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<2; i2++) {
      
      RW *rw = &(orbmat[indxX+i1][indxE+i2]);
      num = rw->num;
      fac = cX[i1]*cE[i2];
      
      if (ngrid != num) {
	cerr << "ResPotOrb::getValues[1]: Oops! ngrid=" 
	     << ngrid << "  num=" << num 
	     << "  indxX=" << indxX+i1 << "/" << numx
	     << "  indxE=" << indxE+i2 << "/" << EX[indxX+i1].size()
	     << endl;
      }
      
      O1  += fac * rw->O1;
      O2  += fac * rw->O2;
    }
  }
  
  // return noboundary;
  return true;
}

bool ResPotOrb::getValues(double rsat, double I1, double I2,
			  double& O1, double& O2,
			  double& Jm, double& dJm,
			  std::complex<double>& Ul, std::complex<double>& dUldE, std::complex<double>& dUldK)
{
  O1 = 0.0;
  O2 = 0.0;
  Jm = 0.0;
  dJm = 0.0;
  Ul = 0.0;
  dUldE = 0.0;
  dUldK = 0.0;
  
  double AMPLITUDE = -4.0*M_PI*MASS*Ylm01(L, M)/rsat;

  int indxX;
  int indxE;
  double cX[2];
  double cE[2];
  
  bool wasok = true;
  
  if (std::isnan(I1) || std::isnan(I2)) {
    cerr << "NaN on input values to getInterp\n";
    wasok = false;
  }
  
  
  bool noboundary;
  getInterp(I1, I2, indxX, indxE, cX, cE, noboundary);
  
  
  for (int i1=0; i1<2; i1++) {
    if (indxE >= (int)EX[indxX+i1].size()-1) {
      cerr << "How did this happen?! [getValue[2]]\n";
      wasok = false;
    }
  }
  
  
  // Compute Ul and derivatives, Jmax, and frequencies
  // -------------------------------------------------
  
  double fac;
  int num;
  
  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<2; i2++) {
      
      RW *rw = &(orbmat[indxX+i1][indxE+i2]);
      num = rw->num;
      fac = cX[i1]*cE[i2];
      
      if (std::isnan(fac)) {
	if (wasok) 
	  cerr << "ResPotOrb::getValues[2]: fac=NaN and was OK!!  cX[" << i1 << "]=" << cX[i1]
	       << "  cE[" << i1 << "][" << i2 << "]=" << cE[i2] 
	       << "  indxX=" << indxX+i1 << "/" << numx
	       << "  indxE=" << indxE+i2 << "/" << EX[i1].size()-1
	       << endl;
	else
	  cerr << "ResPotOrb::getValues[2]: fac=NaN!!  cX[" << i1 << "]=" << cX[i1]
	       << "  cE[" << i1 << "][" << i2 << "]=" << cE[i2] 
	       << "  indxX=" << indxX+i1 << "/" << numx
	       << "  indxE=" << indxE+i2 << "/" << EX[i1].size()-1
	       << endl;
      }
      
      if (ngrid != num) {
	assert( EX[indxX+i1].size() == I1X[indxX+i1].size() );
	cerr << "ResPotOrb::getValues[2]: Oops! ngrid=" 
	     << ngrid << "  num=" << num 
	     << "  i1, i2=" << i1 << ", " << i2
	     << "  indxX=" << indxX+i1 << "/" << numx
	     << "  indxE=" << indxE+i2 << "/" << EX[indxX+i1].size()
	     << endl;
      }
      
      O1  += fac * rw->O1;
      O2  += fac * rw->O2;
      Jm  += fac * rw->Jm;
      dJm += fac * rw->dJm;
      
      // Compute potential transform
      // ---------------------------
      
      W = dWE = dWK = 0.0;
      for (int i=0; i<ngrid; i++) {

	if (rw->r[i]>rsat)
	  W += rw->ff[i]*pow(rw->r[i]/rsat, -(1.0+L));
	else
	  W += rw->ff[i]*pow(rw->r[i]/rsat, L);

	if (rw->r_Ep[i]>rsat)
	  dWE += rw->ff_Ep[i]*pow(rw->r_Ep[i]/rsat, -(1.0+L));
	else
	  dWE += rw->ff_Ep[i]*pow(rw->r_Ep[i]/rsat, L);

	if (rw->r_Em[i]>rsat)
	  dWE -= rw->ff_Em[i]*pow(rw->r_Em[i]/rsat, -(1.0+L));
	else
	  dWE -= rw->ff_Em[i]*pow(rw->r_Em[i]/rsat, L);


	if (rw->r_Kp[i]>rsat)
	  dWK += rw->ff_Kp[i]*pow(rw->r_Kp[i]/rsat, -(1.0+L));
	else
	  dWK += rw->ff_Kp[i]*pow(rw->r_Kp[i]/rsat, L);

	if (rw->r_Km[i]>rsat)
	  dWK -= rw->ff_Km[i]*pow(rw->r_Km[i]/rsat, -(1.0+L));
	else
	  dWK -= rw->ff_Km[i]*pow(rw->r_Km[i]/rsat, L);
      }
  
      Ul +=    fac * W   * AMPLITUDE;
      dUldE += fac * dWE * AMPLITUDE/delE;
      dUldK += fac * dWK * AMPLITUDE/delK;
    }
  }
  
  // return noboundary;
  return true;
}


ResPotOrb::ReturnCode ResPotOrb::coord(double* pos, double* vel,
				 double I1, double I2, double beta,
				 double w1, double w2, double w3)
{
  ofstream out;

  // Linear interpolation coefficients
  // ---------------------------------
  
  int indxX;
  int indxE;
  double cX[2];
  double cE[2];
  
  bool noboundary;
  getInterp(I1, I2, indxX, indxE, cX, cE, noboundary);
  
  for (int i1=0; i1<2; i1++) {
    if (indxE >= (int)EX[indxX+i1].size()-1) {
      cerr << "How did this happen?! [after getInterp]\n";
    }
  }
  
  // Compute radius, f(w1)
  // --------------------
  
  double fac;
  int num;
  
  vector<double> tw(ngrid, 0.0);
  vector<double> tf(ngrid, 0.0);
  vector<double> tr(ngrid, 0.0);
  
  double E = 0.0, test = 0.0;
  double rmin=Rmax, rmax=Rmin;
  
  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<2; i2++) {
      
      RW *rw = &(orbmat[indxX+i1][indxE+i2]);
      num = rw->num;
      fac = cX[i1]*cE[i2];
      test += fac;
      
      if (ngrid != num) {
	cerr << "ResPotOrb::coord[2]: Oops! ngrid=" 
	     << ngrid << "  num=" << num
	     << "  indxX=" << indxX+i1 << "/" << numx
	     << "  indxE=" << indxE+i2 << "/" << EX[indxX+i1].size()
	     << endl;
      }
      
      E += fac * rw->E;
      
      for (int k=0; k<ngrid; k++) {
	tw[k] += fac * rw->w1[k];
	tf[k] += fac * rw->f[k];
	tr[k] += fac * rw->r[k];
      }
      
      rmin = min<double>(rmin, rw->r[0]);
      rmax = max<double>(rmax, rw->r[ngrid-1]);
    }
  }
  
#ifdef DEBUG_DEBUG
  if ( fabs(test-1.0) >= 1.0e-10 ) {
    cout << "Test=" << test << endl;
  }
#endif
  
  assert( fabs(test-1.0)<1.0e-10 );
  
  // Wrap w_1 in [0, 2*pi]
  if (w1>=0.0)
    w1 -=  2.0*M_PI*(int)(0.5*w1/M_PI);
  else
    w1 +=  2.0*M_PI*((int)(-0.5*w1/M_PI) + 1);
  
  double w10 = w1, rv=1.0;
  if (w1>M_PI) {
    w10 = 2.0*M_PI - w1;
    rv = -1.0;
  }
  
  // Compute coordinates
  // -------------------
  
  double r    = odd2(w10, tw, tr);
  double f    = odd2(w10, tw, tf);
  double psi  = w2 - rv*f;
  
  double cosp = cos(psi);
  double sinp = sin(psi);
  
  double cosb = cos(beta);
  double sinb = sin(beta);
  
  double vtot = sqrt(fabs(2.0*(E - halo_model->get_pot(r))));
  double vt   = I2/r;
  
  if (vtot < vt) {		// Adjust r to make orbit valid
    double rcirc;
    ZBrent< double > circ;
    ZBrent< vector<double> > tp;
    ZBrentReturn ret;
    vector<double> param(2);
    param[0] = I2;
    param[1] = E;
    double fmin, fmax;
    const double ftol = 1.0e-3;
    
    if ( (ret=circ.find(get_rc, I2, 0.5*rmin, 2.0*rmax, 1.0e-10, rcirc)) ) {
#ifdef DEBUG_VERBOSE
      cout << "  Rcirc bounds error: val, rmin=" 
	   << get_rc(0.5*rmin, I2) << ", " << rmin
	   << "  val, rmax=" 
	   << get_rc(2.0*rmax, I2) << ", " << rmax << endl;
#endif
      return CoordCirc;
    }
    
    rv = (w1>=M_PI) ? -1.0 : 1.0;
    
    if (w1>1.5*M_PI || w1 < 0.5*M_PI) {	// Below circular orbit radius
      
      fmin = adj_r(Rmin,  param);
      fmax = adj_r(rcirc, param);
      
      if (fmin*fmax >= 0.0) {
	if (fabs(fmin) < ftol*fabs(fmax)) {
	  r = Rmin;
	} else if (fabs(fmax) < ftol*fabs(fmin)) {
	  r = rcirc;
	} else {
#ifdef DEBUG_VERBOSE	
	  cout << "  Inner turning point problem: E=" << E << "  J=" << I2
	       << "  val, r=" << fmin << ", " << Rmin
	       << "  val, rcirc=" << fmax << ", " << rcirc << endl;
#endif
	  return CoordRange;
	}
      } else {			// Above or at circular orbit radius
	
	if ( (ret=tp.find(adj_r, param, Rmin, rcirc, 1.0e-10, r)) ) {
#ifdef DEBUG_VERBOSE	
	  cout << "  Radj inner bounds error: E=" << E << "  J=" << I2
	       << "  val, r=" << adj_r(Rmin, param) 
	       << ", " << Rmin
	       << "  val, rcirc=" << adj_r(rcirc, param) 
	       << ", " << rcirc << endl;
#endif
	  return CoordTP;
	}
	
      }
      
    } else {
      
      fmin = adj_r(rcirc, param);
      fmax = adj_r(Rmax,  param);
      
      if (fmin*fmax >= 0.0) {
	if (fabs(fmin) < ftol*fabs(fmax)) {
	  r = rcirc;
	} else if (fabs(fmax) < ftol*fabs(fmin)) {
	  r = Rmax;
	} else {
#ifdef DEBUG_VERBOSE	
	  cout << "  Outer turning point problem: E=" << E << "  J=" << I2
	       << "  val, r=" << fmin << ", " << Rmin
	       << "  val, rcirc=" << fmax << ", " << rcirc << endl;
#endif
	  return CoordRange;
	}
      } else {
	
	if ( (ret=tp.find(adj_r, param, rcirc, Rmax, 1.0e-10, r)) ) {
#ifdef DEBUG_VERBOSE
	  cout << "  Radj outer bounds error: E=" << E << "  J=" << I2
	       << "  val, rcirc=" << adj_r(rcirc, param) << ", " << rcirc
	       << "  val, r=" << adj_r(Rmax, param) << ", " << Rmax << endl;
#endif
	  return CoordTP;
	}
	
      }
    }
    
    vtot = sqrt(fabs(2.0*(E - halo_model->get_pot(r))));
    vt = I2/r;
    f = 0.0;
    psi = w2;
  }
  
  double vr = 0.0;
  if (fabs(vtot) > fabs(vt)) vr = sqrt(vtot*vtot - vt*vt) * rv;
  
  double cost = sinp*sinb;
  double sint = sqrt(fabs(1.0 - cost*cost));
  
  // Check for excursion beyond bounds
  // ---------------------------------

  if (fabs(cost)>1.0) {
    cost = copysign(1.0, cost);
    sint = 0.0;
  }
  
  double cos3 = cos(w3);
  double sin3 = sin(w3);
  
  // Check for excursion beyond bounds
  // ---------------------------------

  double phi;
  double tmp  = atan2(sin(psi), cos(psi));
  
  if (fabs(sinb)>1.0e-8) {
    
    phi = cosb*cost/(sinb*sint);
    if (fabs(phi)>=1.0) phi = copysign(1.0, phi);
    phi  = asin(phi);
    
    // Phi branch based on Psi
    if (tmp>0.5*M_PI || tmp<-0.5*M_PI) phi = M_PI - phi;
    phi += w3;
    
  } else {
    phi = w3 + tmp;
  }
  
  // Compute positions
  // -----------------
  
  double cosf = cos(phi);
  double sinf = sin(phi);
  
  pos[0] = r * sint*cosf;
  pos[1] = r * sint*sinf;
  pos[2] = r * cost;
  
  // Check final values
  // ------------------

  for (int k=0; k<3; k++)
    if (std::isnan(pos[k]) || std::isinf(pos[k])) return CoordBadPos;


  // Sanity check
  // ------------

  if (std::isnan(vr)  || std::isinf(vr)  || 
      std::isinf(vt)  || std::isinf(vt)  ||
      std::isinf(w3)  || std::isinf(w3)  ||
      std::isinf(psi) || std::isinf(psi)
      )
    {
      pthread_mutex_lock(&iolock);
      out.open(outdir + dbgFILE, ios::app);
      if (out.good()) {
	out <<  "Coord: vel variable is NaN or Inf: vr=" << vr 
	    << " vt=" << vt << " psi=" << psi 
	    << " w3=" << w3 << " r=" << r 
	    << endl;
	out.close();
      }
      pthread_mutex_unlock(&iolock);
    }
  

  // Compute velocities
  // ------------------
  
  double xp = cosp*vr - sinp*vt;
  double yp = sinp*vr + cosp*vt;
  
  vel[0] = xp*cos3 - yp*cosb*sin3;
  vel[1] = xp*sin3 + yp*cosb*cos3;
  vel[2] =           yp*sinb;
  

  // Check final values
  // ------------------

  for (int k=0; k<3; k++)
    if (std::isnan(vel[k]) || std::isinf(vel[k])) return CoordBadVel;
  
  // OK!
  // ---

  return OK;
}


ResPotOrb::ReturnCode ResPotOrb::Update(double dt, 
					vector<double>& Phase, 
					double rsat, double amp,
					double* posI, double* velI,
					double* posO, double* velO)
{
  if (M)
    return Update3(dt, Phase, rsat, amp, posI, velI, posO, velO);
  else
    return Update2(dt, Phase, rsat, amp, posI, velI, posO, velO);
}


ResPotOrb::ReturnCode ResPotOrb::Update2(double dt, 
					 vector<double>& Phase, 
					 double rsat, double amp,
					 double* posI, double* velI,
					 double* posO, double* velO)
{
  ofstream out;
  
  compute_grid();
  
  ReturnCode ret = OK;
  
  if (L1==0 && L2==0) {
    for (int k=0; k<3; k++) {
      posO[k] = posI[k];
      velO[k] = velI[k];
    }
    return UpdateBadL;
  }
  
  // Get action angle coords
  // -----------------------

  double E, W1, W2, W3, F, BETA, PSI;
  double I1, I2, O1, O2;
  ret = coord(posI, velI, E, Kupd, I1, I2, O1, O2, W1, W2, W3, F, BETA, PSI);
  if (ret != OK) return ret;
  
  double betaM, betaP, beta;
  if (BETA-DELTA_B<0.0) {
    betaM = 0.0;
    betaP = DELTA_B;
    beta = 0.5*DELTA_B;
  } else if (BETA+DELTA_B>M_PI) {
    betaM = M_PI - DELTA_B;
    betaP = M_PI;
    beta = M_PI - 0.5*DELTA_B;
  } else {
    betaM = BETA - DELTA_B;
    betaP = BETA + DELTA_B;
    beta = BETA;
  }
  
  std::complex<double> VB = VeeBeta(L, L2, M, beta);
  
  std::complex<double> DVB = 
    (VeeBeta(L, L2, M, betaP) - VeeBeta(L, L2, M, betaM)) / (betaP - betaM);
  
  // Iterative implicit solution
  // ---------------------------
  double Is [4] = {0.0, 0.0, 0.0, 0.0};
  double ws [4] = {0.0, 0.0, 0.0, 0.0};
  double wf [4] = {0.0, 0.0, 0.0, 0.0};
  //                ^    ^    ^    ^
  //                |    |    |    |
  //                |    |    |    \_ Penultimate (for convergence check)
  //                |    |    |
  //                |    |    \_ Latest iteration
  //                |    | 
  //                |    \_ Last step
  //                |
  //                \_ Initial step
  //
  double If, Jm, dJm, dEIs=0.0, dKIs=0.0, dKI1=0.0, dKI2=0.0, I3;
  std::complex<double> Fw, FI, Ul, Ff, dUldE, dUldK, dUldI2, UldVdIs;
  bool done = false;
  
#ifdef DEBUG_DEBUG
  double I10 = I1;
  double I20 = I2;
#endif
  
  if (std::isnan(I1) || std::isnan(I2)) {
    cerr << "Have a cow!\n";
  }
  
  // Transformation to slow-fast variables
  // -------------------------------------
  
  ws[2]  = ws[0]  = W1*L1 + W2*L2 + (W3 - Phase[0])*M;
  if (L2) {
    Is[2] = Is[0] = I2/L2;
    wf[2] = wf[0] = W1;
    If = I1 - Is[0]*L1;
  } else {
    Is[2] = Is[0] = I1/L1;
    wf[2] = wf[0] = W2;
    If = I2 - Is[0]*L2;
  }
  
  I3 = I2*cos(BETA);
  
  double Omega = (Phase[2] - Phase[0])/dt;

  int i;
  for (i=0; i<ITMAX; i++) {
    
    // For convergence test
    // --------------------

    ws[3] = ws[1]; 
    wf[3] = wf[1];
    Is[3] = Is[1];
    
    
    // Save previous step
    // ------------------

    ws[1] = ws[2];
    wf[1] = wf[2];
    Is[1] = Is[2];
    
    // For force evaluation
    // --------------------

    if (second_order) {
      // Second order
      Is[2] = 0.5*(Is[2] + Is[0]);
      ws[2] = 0.5*(ws[2] + ws[0]);
      wf[2] = 0.5*(wf[2] + wf[0]);
      
    } else {
      // First order
      Is[2] = Is[2];
      ws[2] = ws[2];
      wf[2] = wf[2];
    }
    
    // Canonical transform
    // -------------------

    if (L2) {
      I1 = If + Is[2]*L1;
      I2 = Is[2]*L2;
    } else {
      I1 = Is[2]*L1;
      I2 = If + Is[2]*L2;
    }
    
    getValues(rsat, I1, I2, O1, O2, Jm, dJm, Ul, dUldE, dUldK);

    // Sanity check
    // ------------

    if (std::isnan(I1) || std::isnan(I2)) {
      cerr << "I1 or I2 is NaN: Is0=" 
	   << Is[0] << " Is=" << Is << " If=" 
	   << If << " Is_2=" << Is[2] << " i=" 
	   << i << endl;
      
      pthread_mutex_lock(&iolock);
      out.open(outdir+dbgFILE, ios::app);
      if (out.good()) {
	out <<  "I1 or I2 is NaN: Is0=" 
	    << Is[0] << " Is1=" << Is[1] << " If=" 
	    << If << " Is2=" << Is[2] << " i=" 
	    << i << endl;
	out.close();
      }
      pthread_mutex_unlock(&iolock);
    }
    
    
    dUldI2 = Ul * DVB * amp / (tan(beta)*I2);
    UldVdIs = dUldI2 * static_cast<double>(L2);
    Ul *= VB * amp;
    dUldE *= VB * amp;
    dUldK *= VB * amp;
    
    dEIs = O1*L1 + O2*L2;
    dKIs = 1.0/Jm*L2;
    dKI1 = -I2*dJm/(Jm*Jm)*O1;
    dKI2 = 1.0/Jm - I2*dJm/(Jm*Jm)*O2;
    
    Fw = O1*L1 + O2*L2 - Omega*M +
      (dUldE*dEIs + dUldK*dKIs + UldVdIs)*exp(I*ws[2]);
    FI = -I*Ul*exp(I*ws[2]);
    if (L2)
      Ff = O1 + (dUldE*O1 + dUldK*dKI1)*exp(I*ws[2]);
    else
      Ff = O2 + (dUldE*O2 + dUldK*dKI2 + dUldI2)*exp(I*ws[2]);
    
    // Sanity check
    // ------------

    if (
	std::isnan(Fw.real()) || std::isnan(FI.real()) ||
	std::isnan(Ff.real())
	) {
      cerr << "Fw or FI is NaN, dJm=" << dJm 
	   << " Ul="	<< Ul 
	   << " dUldE="	<< dUldE 
	   << " dUldK="	<< dUldK 
	   << " dEIs="	<< dEIs 
	   << " dKIs="	<< dKIs 
	   << " O1="	<< O1 
	   << " O2="	<< O2 
	   << " Omega="	<< Omega 
	   << " P0="	<< Phase[0]
	   << " P2="	<< Phase[2]
	   << " dt="	<< dt
	   << " ws="	<< ws[1]
	   << " ws0="	<< ws[0]
	   << " ws2="	<< ws[2]
	   << " i="	<< i << endl;
      
      pthread_mutex_lock(&iolock);
      out.open(outdir+dbgFILE, ios::app);
      if (out.good()) {
	out  << "Fw or FI is NaN, dJm=" << dJm 
	     << " Ul="	<< Ul 
	     << " dUldE="	<< dUldE 
	     << " dUldK="	<< dUldK 
	     << " dEIs="	<< dEIs 
	     << " dKIs="	<< dKIs 
	     << " O1="	<< O1 
	     << " O2="	<< O2 
	     << " Omega="	<< Omega 
	     << " P0="	<< Phase[0]
	     << " P2="	<< Phase[2]
	     << " dt="	<< dt
	     << " ws="	<< ws[1]
	     << " ws0="	<< ws[0]
	     << " ws2="	<< ws[2]
	     << " i="	<< i << endl;
	out.close();
      }
      pthread_mutex_unlock(&iolock);
    }
    
    // Update
    // ------

    ws[2] = ws[0] + dt*Fw.real();
    Is[2] = Is[0] + dt*FI.real();
    wf[2] = wf[0] + dt*Ff.real();
    
    if (std::isnan(ws[2])) {
      cerr << "ws2 is NaN, Fw=" << Fw.real()
	   << " ws0=" << ws[0]
	   << " dt=" << dt
	   << " i="	<< i << endl;
      
      pthread_mutex_lock(&iolock);
      out.open(outdir+dbgFILE, ios::app);
      if (out.good()) {
	out  << "ws2 is NaN, Fw=" << Fw.real()
	     << " ws0=" << ws[0]
	     << " dt=" << dt
	     << " i="	<< i << endl;
	out.close();
      }
      pthread_mutex_unlock(&iolock);
    }
    
    
    // Check for convergence
    // ---------------------

    if (fabs(ws[1]-ws[2])<TOLITR*dt && 
	fabs(Is[1]-Is[2])/(fabs(Is[2])+1.0e-10)<TOLITR*dt &&
	fabs(wf[1]-wf[2])<TOLITR*dt
	) done = true;
    
    // Limit 2-cycle detection
    // -----------------------

    if (i>3 &&
	fabs(ws[3]-ws[2])<TOLITR*dt && 
	fabs(Is[3]-Is[2])/(fabs(Is[2])+1.0e-10)<TOLITR*dt) done = true;
    
    if (done) break;
    
  }
  
  if (!done) {
    pthread_mutex_lock(&iolock);
    out.open(outdir+dbgFILE, ios::app);
    if (out.good()) {
      out << "Update iteration: "
	  << "Phase, E, K, I1, I2, DI, Dw, Ul, dUldE, dUldK, dEIs, dKIs = " 
	  << Phase[1]
	  << ", " << E
	  << ", " << Kupd
	  << ", " << I1
	  << ", " << I2
	  << ", " << Is[2]-Is[1]
	  << ", " << ws[2]-ws[1]
	  << ", " << Is[2]-Is[3]
	  << ", " << ws[2]-ws[3]
	  << ", " << Ul
	  << ", " << dUldE
	  << ", " << dUldK
	  << ", " << dEIs
	  << ", " << dKIs
	  << endl;
      out.close();
    }
    pthread_mutex_unlock(&iolock);
    ret = UpdateIterate;
  }
  
  // Canonical transformation from slow-fast to action-angle
  // -------------------------------------------------------

  if (L2) {
    W1 = wf[2];
    W2 = (ws[2] - W1*L1 + Phase[2]*M)/L2;
    I1 = Is[2]*L1 + If;
    I2 = Is[2]*L2;
  } else {
    W2 = wf[2];
    W1 = (ws[2] - W2*L2 + Phase[2]*M)/L2;
    I1 = Is[2]*L1;
    I2 = Is[2]*L2 + If;
  }
  
  double cosb = I3/I2;
  cosb = min<double>(cosb,  1.0);
  cosb = max<double>(cosb, -1.0);
  BETA = acos(cosb);
  
#ifdef DEBUG_DEBUG
  if (fabs(I10-I1)>1.0e-3*fabs(I10) || fabs(I20-I2)>1.0e-3*fabs(I20)) {
    cout << setw(15) << I10
	 << setw(15) << I1
	 << setw(15) << I20
	 << setw(15) << I2
	 << setw(15) << fabs(I10 - I1)
	 << setw(15) << fabs(I20 - I2)
	 << endl;
  }
#endif
  
  // Get new Cartesion phase space
  // -----------------------------

  ret = coord(posO, velO, I1, I2, BETA, W1, W2, W3); 
  if (ret != OK) {
    for (int k=0; k<3; k++) {
      posO[k] = posI[k];
      velO[k] = velI[k];
    }
  }
  
  // Return last completion code
  // ---------------------------

  return ret;
}


ResPotOrb::ReturnCode ResPotOrb::Update3(double dt, 
					 vector<double>& Phase, 
					 double rsat, double amp, 
					 double* posI, double* velI,
					 double* posO, double* velO)
{
  std::ofstream out;
  
  compute_grid();
  
  ReturnCode ret;
  
  if (L1==0 && L2==0) {
    for (int k=0; k<3; k++) {
      posO[k] = posI[k];
      velO[k] = velI[k];
    }
    return UpdateBadL;
  }
  
  // Get action angle coords
  // -----------------------
  double E, W1, W2, W3, F, BETA, PSI;
  double I1, I2, O1, O2;
  if ((ret=coord(posI, velI, E, Kupd, I1, I2, O1, O2, W1, W2, W3, F, BETA, PSI)) != OK) return ret;
  
  double betaM, betaP, beta;
  if (BETA-DELTA_B<0.0) {
    betaM = 0.0;
    betaP = DELTA_B;
    beta = 0.5*DELTA_B;
  } else if (BETA+DELTA_B>M_PI) {
    betaM = M_PI - DELTA_B;
    betaP = M_PI;
    beta = M_PI - 0.5*DELTA_B;
  } else {
    betaM = BETA - DELTA_B;
    betaP = BETA + DELTA_B;
    beta = BETA;
  }
  
  std::complex<double> VB = VeeBeta(L, L2, M, beta);
  
  std::complex<double> DVB = 
    (VeeBeta(L, L2, M, betaP) - VeeBeta(L, L2, M, betaM)) / (betaP - betaM);
  
  // Iterative implicit solution
  // ---------------------------
  double Is [4] = {0.0, 0.0, 0.0, 0.0};
  double ws [4] = {0.0, 0.0, 0.0, 0.0};
  double wf1[4] = {0.0, 0.0, 0.0, 0.0};
  double wf2[4] = {0.0, 0.0, 0.0, 0.0};
  //                ^    ^    ^    ^
  //                |    |    |    |
  //                |    |    |    \_ Penultimate (for convergence check)
  //                |    |    |
  //                |    |    \_ Latest iteration
  //                |    | 
  //                |    \_ Last step
  //                |
  //                \_ Initial step
  //
  double If1, If2;
  double Jm, dJm, dEIs=0.0, dKIs=0.0, dKI1=0.0, dKI2=0.0, I3, I30;
  std::complex<double> Fw, FI, Ul, F1, F2, dUldE, dUldK, dUldb;
  bool done = false;
  
#ifdef DEBUG_DEBUG
  double I10 = I1;
  double I20 = I2;
#endif
  
  if (std::isnan(I1) || std::isnan(I2)) {
    cerr << "Have a cow!\n";
  }
  
  // Transformation to slow-fast variables
  // -------------------------------------
  
  ws[2]  = ws[0]  = W1*L1 + W2*L2 + (W3 - Phase[0])*M;
  wf1[2] = wf1[0] = W1;
  wf2[2] = wf2[0] = W2;
  
  I30 = I2*cos(beta);
  
  Is[2] = Is[0] = I30/M;
  If1 = I1 - Is[0]*L1;
  If2 = I2 - Is[0]*L2;
  
  
  double Omega = (Phase[2] - Phase[0])/dt;

  
  int i;
  for (i=0; i<ITMAX; i++) {
    
    // For convergence test
    // --------------------
    ws [3] = ws [1]; 
    wf1[3] = wf1[1];
    wf2[3] = wf2[1];
    Is [3] = Is [1];
    
    
    // Save previous step
    // ------------------
    ws [1] = ws [2];
    wf1[1] = wf1[2];
    wf2[1] = wf2[2];
    Is [1] = Is [2];
    
    // For force evaluation
    // --------------------
    if (second_order) {
      // Second order
      Is[2]  = 0.5*(Is[2]  + Is[0]);
      ws[2]  = 0.5*(ws[2]  + ws[0]);
      wf1[2] = 0.5*(wf1[2] + wf1[0]);
      wf2[2] = 0.5*(wf1[2] + wf1[0]);
      
    } else {
      // First order
      Is[2] = Is[2];
      ws[2] = ws[2];
      wf1[2] = wf1[2];
      wf2[2] = wf2[2];
    }
    
    // Canonical transform
    // -------------------
    I1 = If1 + Is[2]*L1;
    I2 = If2 + Is[2]*L2;
    
    getValues(rsat, I1, I2, O1, O2, Jm, dJm, Ul, dUldE, dUldK);

    // Sanity check
    // ------------
    if (std::isnan(I1) || std::isnan(I2)) {
      cerr << "I1 or I2 is NaN: Is0=" 
	   << Is[0] << " Is=" << Is << " If1=" 
	   << If1 << " If2=" << If2 << " Is_2=" << Is[2] << " i=" 
	   << i << endl;
      
      pthread_mutex_lock(&iolock);
      out.open(outdir+dbgFILE, ios::app);
      if (out.good()) {
	out <<  "I1 or I2 is NaN: Is0=" 
	    << Is[0] << " Is1=" << Is[1] << " If1=" 
	    << If1 << " If2=" << If2 << " Is2=" << Is[2] << " i=" 
	    << i << endl;
	out.close();
      }
      pthread_mutex_unlock(&iolock);
    }
    
    
    dUldb = -Ul * DVB * amp * static_cast<double>(M) / (sin(beta)*I2);
    Ul *= VB * amp;
    dUldE *= VB * amp;
    dUldK *= VB * amp;
    
    dEIs = O1*L1 + O2*L2;
    dKIs = 1.0/Jm*L2;
    dKI1 = -I2*dJm/(Jm*Jm)*O1;
    dKI2 = 1.0/Jm - I2*dJm/(Jm*Jm)*O2;
    
    Fw = O1*L1 + O2*L2 - Omega*M + dUldb*exp(I*ws[2]);
    FI = -I*Ul*exp(I*ws[2]);
    F1 = O1 + (dUldE*O1 + dUldK*dKI1)*exp(I*ws[2]);
    F2 = O2 + (dUldE*O2 + dUldK*dKI2)*exp(I*ws[2]);
    
    // Sanity check
    // ------------
    if (
	std::isnan(Fw.real()) || std::isnan(FI.real()) ||
	std::isnan(F1.real()) || std::isnan(F2.real())
	) {
      cerr << "Fw or FI is NaN, dJm=" << dJm 
	   << " Ul="	<< Ul 
	   << " dUldE="	<< dUldE 
	   << " dUldK="	<< dUldK 
	   << " dEIs="	<< dEIs 
	   << " dKIs="	<< dKIs 
	   << " O1="	<< O1 
	   << " O2="	<< O2 
	   << " Omega="	<< Omega 
	   << " P0="	<< Phase[0]
	   << " P0="	<< Phase[2]
	   << " dt="	<< dt
	   << " ws="	<< ws[1]
	   << " ws0="	<< ws[0]
	   << " ws2="	<< ws[2]
	   << " i="	<< i << endl;
      
      pthread_mutex_lock(&iolock);
      out.open(outdir+dbgFILE, ios::app);
      if (out.good()) {
	out  << "Fw or FI is NaN, dJm=" << dJm 
	     << " Ul="	<< Ul 
	     << " dUldE="	<< dUldE 
	     << " dUldK="	<< dUldK 
	     << " dEIs="	<< dEIs 
	     << " dKIs="	<< dKIs 
	     << " O1="	<< O1 
	     << " O2="	<< O2 
	     << " Omega="	<< Omega 
	     << " P0="	<< Phase[0]
	     << " P0="	<< Phase[2]
	     << " dt="	<< dt
	     << " ws="	<< ws[1]
	     << " ws0="	<< ws[0]
	     << " ws2="	<< ws[2]
	     << " i="	<< i << endl;
	out.close();
      }
      pthread_mutex_unlock(&iolock);
    }
    
    // Update
    // ------
    ws[2]  = ws[0]  + dt*Fw.real();
    Is[2]  = Is[0]  + dt*FI.real();
    wf1[2] = wf1[0] + dt*F1.real();
    wf2[2] = wf2[0] + dt*F2.real();
    
    if (std::isnan(ws[2])) {
      cerr << "ws2 is NaN, Fw=" << Fw.real()
	   << " ws0=" << ws[0]
	   << " dt=" << dt
	   << " i="	<< i << endl;
      
      pthread_mutex_lock(&iolock);
      out.open(outdir+dbgFILE, ios::app);
      if (out.good()) {
	out  << "ws2 is NaN, Fw=" << Fw.real()
	     << " ws0=" << ws[0]
	     << " dt=" << dt
	     << " i="	<< i << endl;
	out.close();
      }
      pthread_mutex_unlock(&iolock);
    }
    
    
    // Check for convergence
    // ---------------------
    if (fabs(ws[1]-ws[2])<TOLITR*dt && 
	fabs(Is[1]-Is[2])/(fabs(Is[2])+1.0e-10)<TOLITR*dt &&
	fabs(wf1[1]-wf1[2])<TOLITR*dt && 
	fabs(wf2[1]-wf2[2])<TOLITR*dt
	) done = true;
    
    // Limit 2-cycle detection
    // -----------------------
    if (i>3 &&
	fabs(ws[3]-ws[2])<TOLITR*dt && 
	fabs(Is[3]-Is[2])/(fabs(Is[2])+1.0e-10)<TOLITR*dt) done = true;
    
    if (done) break;
    
  }
  
  if (!done) {
    pthread_mutex_lock(&iolock);
    out.open(outdir+dbgFILE, ios::app);
    if (out.good()) {
      out << "Update iteration: "
	  << "Phase, E, K, I1, I2, DI, Dw, Ul, dUldE, dUldK, dEIs, dKIs = " 
	  << Phase[1]
	  << ", " << E
	  << ", " << Kupd
	  << ", " << I1
	  << ", " << I2
	  << ", " << Is[2]-Is[1]
	  << ", " << ws[2]-ws[1]
	  << ", " << Is[2]-Is[3]
	  << ", " << ws[2]-ws[3]
	  << ", " << Ul
	  << ", " << dUldE
	  << ", " << dUldK
	  << ", " << dEIs
	  << ", " << dKIs
	  << endl;
      out.close();
    }
    pthread_mutex_unlock(&iolock);
    ret = UpdateIterate;
  }
  
  // Canonical transformation from slow-fast to action-angle
  // -------------------------------------------------------
  W1 = wf1[2];
  W2 = wf2[2];
  W3 = (ws[2] - W1*L1 - W2*L2 + Phase[2]*M)/M;
  I1 = Is[2]*L1 + If1;
  I2 = Is[2]*L2 + If2;
  I3 = Is[2]*M;
  
  double cosb = I3/I2;
  cosb = min<double>(cosb,  1.0);
  cosb = max<double>(cosb, -1.0);
  BETA = acos(cosb);
  
#ifdef DEBUG_DEBUG
  if (fabs(I10-I1)>1.0e-3*fabs(I10) || fabs(I20-I2)>1.0e-3*fabs(I20)) {
    cout << setw(15) << I10
	 << setw(15) << I1
	 << setw(15) << I20
	 << setw(15) << I2
	 << setw(15) << fabs(I10 - I1)
	 << setw(15) << fabs(I20 - I2)
	 << endl;
  }
#endif
  
  // Get new Cartesion phase space
  // -----------------------------
  if ((ret=coord(posO, velO, I1, I2, BETA, W1, W2, W3)) != OK) {
    for (int k=0; k<3; k++) {
      posO[k] = posI[k];
      velO[k] = velI[k];
    }
  }
  
  // Return last completion code
  // ---------------------------
  return ret;
}


ResPotOrb::ReturnCode ResPotOrb::Force(double dt, vector<double>& Phase, 
				       double rsat, double amp,
				       double* pos, double* vel, double* acc)
{
  ReturnCode ret = OK;
  double pos0[3], vel0[3], pos2[3], vel2[3];
  double E, K, W1, W2, W3, F, BETA, PSI, I1, I2, O1, O2;
  
  // Get action angle coords
  ret = coord(pos, vel, E, K, I1, I2, O1, O2, W1, W2, W3, F, BETA, PSI);
  
  if (ret != OK) {
    for (int k=0; k<3; k++) acc[k] = 0.0;
    return ret;
  }
  
  // Get phase space update without perturbation
  // -------------------------------------------
  W1 += O1*dt;
  W2 += O2*dt;
  ret = coord(pos0, vel0, I1, I2, BETA, W1, W2, W3);
  
  if (ret != OK) {
    for (int k=0; k<3; k++) acc[k] = 0.0;
    return ret;
  }
  
  // Get phase space update with perturbation
  // ----------------------------------------
  ret = Update(dt, Phase, rsat, amp, pos, vel, pos2, vel2);
  
  if (ret != OK) {
    for (int k=0; k<3; k++) acc[k] = 0.0;
    return ret;
  }
  
  // Effective acceleration
  // ----------------------
  for (int k=0; k<3; k++) acc[k] = (vel2[k] - vel0[k])/dt;
  
  return ret;
}


void ResPotOrb::check_rw(RW& rw)
{
  for (unsigned i=0; i<rw.r.size(); i++)
    if (std::isnan(rw.r[i])) cout << "RW error: r nan, i=" << i << endl;
  
  for (unsigned i=0; i<rw.w1.size(); i++)
    if (std::isnan(rw.w1[i])) cout << "RW error: w1 nan, i=" << i << endl;
  
  for (unsigned i=0; i<rw.f.size(); i++)
    if (std::isnan(rw.f[i])) cout << "RW error: f nan, i=" << i << endl;
  
  if (std::isnan(rw.O1))	cout << "RW error: O1 nan" << endl;
  if (std::isnan(rw.O2))	cout << "RW error: O2 nan" << endl;
  if (std::isnan(rw.Jm))	cout << "RW error: Jm nan" << endl;
  if (std::isnan(rw.dJm))	cout << "RW error: dJm nan" << endl;
  if (std::isnan(rw.E))	        cout << "RW error: E nan" << endl;
  if (std::isnan(rw.K))	        cout << "RW error: K nan" << endl;
  if (std::isnan(rw.I1))	cout << "RW error: I1 nan" << endl;
}

