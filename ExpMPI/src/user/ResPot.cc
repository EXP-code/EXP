#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <values.h>
#include <ResPot.H>
#include <ZBrent.H>
#include <localmpi.h>

#include <pthread.h>  
pthread_mutex_t iolock = PTHREAD_MUTEX_INITIALIZER;

double ResPot::ALPHA = 0.25;
double ResPot::DELTA = 0.01;
double ResPot::DELE = 0.001;
double ResPot::DELK = 0.001;
double ResPot::TOLITR = 1.0e-8;
int ResPot::NUME = 200;
int ResPot::NUMX = 200;
int ResPot::RECS = 100;
int ResPot::ITMAX = 50;
KComplex ResPot::I(0.0, 1.0);


ResPot::ResPot(AxiSymModel *mod, AxiSymBiorth *bio, 
	       int l, int m, int l1, int l2, int nmax)
{
  halo_model = mod;
  halo_ortho = bio;
  
  grid_computed = false;
  second_order = true;
  
  L = l;
  M = m;
  L1 = l1;
  L2 = l2;
  NMAX = nmax;
  
  Rmin = halo_model->get_min_radius();
  Rmax = halo_model->get_max_radius();
  Emax = halo_model->get_pot(Rmax)*(1.0+DELTA);
  Emin = halo_model->get_pot(Rmin)*(1.0-DELTA);
  Kmin = DELTA;
  Kmax = 1.0-DELTA;
  
  // SphericalOrbit::RMAXF=1.0;
  orb = new SphericalOrbit(halo_model);
  orb->set_numerical_params(RECS);
  orb->set_biorth(*halo_ortho, L, NMAX);
  
  compute_grid();
}

ResPot::~ResPot()
{
  delete orb;
}

AxiSymModel *hmod;

// Circular orbit
//
static double get_rc(double r, double J)
{
  return J*J - hmod->get_mass(r)*r;
}

// Turning points
//
static double adj_r(double r, vector<double> p)
{
  return p[1] - 0.5*p[0]*p[0]/(r*r) - hmod->get_pot(r);
}


void ResPot::compute_grid() 
{
  if (grid_computed) return;
  
  double E, K;
  Vector t(1, NMAX);
  
  // Ang mom of minimum radius circular orbit
  minJ = Rmin*sqrt(halo_model->get_dpot(Rmin));
  
  // Ang mom of maximum radius circular orbit
  maxJ = Rmax*sqrt(halo_model->get_dpot(Rmax));
  
  // For derivatives
  double delE = (Emax - Emin)*0.5*DELE;
  double delK = (Kmax - Kmin)*0.5*DELK;
  
  double J, rcirc, Ecirc;
  int ibeg;
  
  dX = 1.0/(NUMX-1);
  double dE = (Emax - Emin)/(NUME-1);
  
  hmod = halo_model;
  
  int iret;
  ZBrent<double> circ;
  
  for (int s=0; s<NUMX; s++) {
    
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
    
    ibeg = (int)( (Ecirc - Emin)/dE );
    
    for (int i=ibeg; i<NUME; i++) {
      
      if (i==ibeg)
	E = Ecirc;
      else
	E = Emin + dE*i;
      
      orb->new_orbit(E, 0.5);
      K = J/orb->Jmax();
      K = max<double>(K, Kmin);
      K = min<double>(K, Kmax);
      
      RW rw;
      
      orb->new_orbit(E, K);
      orb->set_biorth(*halo_ortho, L, NMAX);
      orb->pot_trans(L1, L2, t);
      for (int n=0; n<NMAX; n++) rw.W.push_back(t[n+1]);
      
      struct ANGLE_GRID * grid = orb->get_angle_grid();
      
      for (int j=0; j<grid->num; j++) {
	rw.r.push_back(grid->r[1][j]);
	rw.w1.push_back(grid->w1[1][j]);
	rw.f.push_back(grid->f[1][j]);
      }
      rw.num = grid->num;
      rw.O1 = orb->get_freq(1);
      rw.O2 = orb->get_freq(2);
      rw.I1 = orb->get_action(1);
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
      orb->set_biorth(*halo_ortho, L, NMAX);
      orb->pot_trans(L1, L2, t);
      for (int n=0; n<NMAX; n++) rw.dWE.push_back(t[n+1]);
      rw.dJm = orb->Jmax();
      
      orb->new_orbit(Em, K);
      orb->set_biorth(*halo_ortho, L, NMAX);
      orb->pot_trans(L1, L2, t);
      for (int n=0; n<NMAX; n++) rw.dWE[n] = 0.5*(rw.dWE[n] - t[n+1])/delE;
      rw.dJm = 0.5*(rw.dJm - orb->Jmax())/delE;
      
      // Kappa deriv
      orb->new_orbit(E, Kp);
      orb->set_biorth(*halo_ortho, L, NMAX);
      orb->pot_trans(L1, L2, t);
      for (int n=0; n<NMAX; n++) rw.dWK.push_back(t[n+1]);
      
      orb->new_orbit(E, Km);
      orb->set_biorth(*halo_ortho, L, NMAX);
      orb->pot_trans(L1, L2, t);
      for (int n=0; n<NMAX; n++) rw.dWK[n] = 0.5*(rw.dWK[n] - t[n+1])/delK;
      
      // Finally, store the data for this phase space point
      orbvec.push_back(rw);
      check_rw(rw);
    }
    
    if (EJ1.size()<2) break;
    
    EX.push_back(EJ1);
    I1X.push_back(IJ1);
    
    orbmat.push_back(orbvec);
  }
  
  numx = EX.size();
  ngrid = orb->get_angle_grid()->num;
  
  grid_computed = true;
}


bool ResPot::coord(double* ps1, double* vel,
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
  
  if (r>halo_model->get_max_radius()) return false;
  
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
  
  int indxE[2];
  double cE[2][2];
  for (int i1=0; i1<2; i1++) {
    indxE[i1] = Vlocate(E, EX[indxX+i1]);
    indxE[i1] = min<int>(indxE[i1], EX[indxX+i1].size()-2);
    indxE[i1] = max<int>(indxE[i1], 0);
    
    cE[i1][0] = (EX[indxX+i1][indxE[i1]+1] - E)/
      (EX[indxX+i1][indxE[i1]+1] - EX[indxX+i1][indxE[i1]]);
    
    // Bounds: keep it on the grid else unphsical!
    cE[i1][0] = max<double>(cE[i1][0], 0.0);
    cE[i1][0] = min<double>(cE[i1][0], 1.0);
    
    cE[i1][1] = 1.0 - cE[i1][0];
    
  }
  
  
  // Compute angles
  // --------------
  
  double fac;
  int num;
  
  vector<double> tw(ngrid, 0.0);
  vector<double> tf(ngrid, 0.0);
  vector<double> tr(ngrid, 0.0);
  
  I1 = K = O1 = O2 = 0.0;
  
  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<2; i2++) {
      
      RW *rw = &(orbmat[indxX+i1][indxE[i1]+i2]);
      num = rw->num;
      fac = cX[i1]*cE[i1][i2];
      
      if (ngrid != num) {
	cerr << "ResPot::coord[1]: Oops! ngrid=" 
	     << ngrid << "  num=" << num 
	     << "  indxX=" << indxX+i1 << "/" << numx
	     << "  indxE=" << indxE[i1]+i2 << "/" << EX[indxX+i1].size()
	     << endl;
      }
      
      I1 += fac * rw->I1;
      K  += fac * rw->K;
      O1 += fac * rw->O1;
      O2 += fac * rw->O2;
      
      for (int k=0; k<ngrid; k++) {
	tw[k] += fac * rw->w1[k];
	tf[k] += fac * rw->f[k];
	tr[k] += fac * rw->r[k];
      }
      
    }
  }
  
  double w1 = odd2(r, tr, tw);
  double f  = odd2(r, tr, tf);
  
  if (rv < 0.0) {
    w1 = 2.0*M_PI - w1;
    f *= -1.0;
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
  
  return true;
}


//
// Coordindate mapping to unit interval
//

double ResPot::xJ(double J, double Jmin, double Jmax)
{
  if (J>=Jmax) return 1.0;
  if (J<=Jmin) return 0.0;
  return pow( (J-Jmin)/(Jmax-Jmin), ALPHA );
}

double ResPot::Jx(double x, double Jmin, double Jmax)
{
  if (x>=1.0) return Jmax;
  if (x<=0.0) return Jmin;
  return Jmin + (Jmax - Jmin)*pow(x, 1.0/ALPHA);
}


double ResPot::dxJ(double J, double Jmin, double Jmax)
{
  if (J>=Jmax) return 1.0/(Jmax-Jmin);
  if (J<=0.0) return MAXFLOAT;
  return pow( (J-Jmin)/(Jmax-Jmin), ALPHA-1.0 )/(Jmax - Jmin);
}


void ResPot::getInterp(double I1, double I2, int& indxX, int indxE[2],
		       double cX[2], double cE[2][2], bool& noboundary)
{
  // Linear interpolation coefficients
  // ---------------------------------
  
  double X = xJ(I2, minJ, maxJ);
  
  indxX = (int)( X/dX );
  indxX = min<int>(indxX, numx-2);
  indxX = max<int>(indxX, 0);
  
  cX[0] = (dX*(indxX+1) - X)/dX;
  cX[1] = 1.0 - cX[0];
  
  if (isnan(cX[0]) || isnan(cX[1])) {
    cerr << "getInterp: cX is NaN, X=" << X 
	 << "  I1=" << I1
	 << "  I2=" << I2 
	 << "  minJ=" << minJ
	 << "  maxJ=" << maxJ
	 << endl;
  }

  int numE[2];
  
  // Upper bound
  for (int i1=0; i1<2; i1++) {
    numE[i1]  = I1X[indxX+i1].size()-1;
    indxE[i1] = Vlocate(I1, I1X[indxX+i1]) + 1;
    indxE[i1] = min<int>(indxE[i1], numE[i1]);
    indxE[i1] = max<int>(indxE[i1], 0);
  }
  
  // Get offset from outer energy boundary
  
  int topoff = min<int>(numE[0] - indxE[0], numE[1] - indxE[1]);
  assert( topoff >= 0);
  
  // Search down through energy grid for match
  
  bool ok = false;
  // Simple case: evenly spaced energy grid
  // --------------------------------------
  int minE = min<int>(numE[0], numE[1]) - 1;
  
  double I1cur=0.0;
  for (int i1=0; i1<2; i1++) I1cur += cX[i1]*I1X[indxX+i1][numE[i1]-topoff];
  double I1last = I1cur;

  // Over top of grid
  // ----------------
  if (I1>I1cur) {
    assert( topoff==0 );
    for (int i1=0; i1<2; i1++) {
      indxE[i1] = I1X[indxX+i1].size()-2;
      cE[i1][0] = 0.0;
      cE[i1][1] = 1.0;
    }
    
    noboundary = false;

    return;
  }
  
  // At bottom of grid
  // -----------------
  if (topoff>minE) {
    for (int i1=0; i1<2; i1++) {
      indxE[i1] = 0;
      cE[i1][0] = 1.0;
      cE[i1][1] = 0.0;
    }
    
    noboundary = false;

    return;
  }

  for (int i=topoff+1; i<=minE; i++) {
    assert( i>0 );
    I1last = I1cur;
    I1cur = 0.0;
    for (int i1=0; i1<2; i1++) {
      assert( numE[i1]-i >= 0 );
      I1cur += cX[i1]*I1X[indxX+i1][numE[i1]-i];
    }
    if ( (I1-I1cur)*(I1-I1last) <= 0.0 ) {
      cE[1][0] = cE[0][0] = (I1last - I1)/(I1last - I1cur);
      cE[1][1] = cE[0][1] = 1.0 - cE[0][0];
      indxE[0] = numE[0] - i;
      indxE[1] = numE[1] - i;
      ok = true;
      break;
    }
  }
  
  noboundary = true;

  double I1tst1 = I1last;
  double I1tst2 = I1cur;

  if (!ok) {
    
    noboundary = false;

    assert( numE[0]-minE-1 >= 0 );

    assert( numE[1]-minE-1 == 0 );

    assert(
	   (EX[indxX+1][0] >= EX[indxX][numE[0]-minE-1]) &&
	   (EX[indxX+1][0] <= EX[indxX][numE[0]-minE+0])
	   );
    
    // Special case #1: I1 for energies larger
    //      than the max-min energy endpoint
   
    I1last = I1cur;
    I1cur = 
      cX[0]*odd2(EX[indxX+1][0], EX[indxX], I1X[indxX]) +
      cX[1]*I1X[indxX+1][0];
    
    if ( (I1-I1cur)*(I1-I1last) <= 0.0 ) {

      indxE[0] = numE[0]-minE-1;

      if (fabs(I1last - I1X[indxX][indxE[0]]) < 1.0e-8) {
	cerr << "Interp denom=0 in Special Case #1: I1last=" << I1last
	     << "  I1cur=" << I1cur
	     << "  I1tst1=" << I1tst1
	     << "  I1tst2=" << I1tst2
	     << "  I1prv=" << I1X[indxX][numE[0]-minE]
	     << "  I1min=" << I1X[indxX][numE[0]-minE-1]
	     << "  I1pen=" << I1X[indxX+1][numE[0]-minE]
	     << "  I1int=" << odd2(EX[indxX+1][0], EX[indxX], I1X[indxX])
	     << "  I1=" << I1
	     << "  dI1top=" << I1-I1last
	     << "  dI1bot=" << I1cur-I1last
	     << "  indxE=" << indxE[0]
	     << "  minE=" << minE
	     << "  numE{0]=" << numE[0]
	     << "  numE[1]=" << numE[1]
	     << "  topoff=" << topoff
	     << "  loc indx=" << Vlocate(EX[indxX+1][0], EX[indxX])
	     << "  cX[0]=" << cX[0] << "  CX[1]=" << cX[1]
	     << endl;
      }

      if (fabs(I1last - I1X[indxX][indxE[0]]) > 1.0e-8)
	cE[0][0] = (I1last - I1)/(I1last - I1X[indxX][indxE[0]]);
      else
	cE[0][0] = (EX[indxX][indxE[0]+1] - EX[indxX+1][0])/(EX[indxX][indxE[0]+1] - EX[indxX][indxE[0]]);
      cE[0][0] = min<double>(cE[0][0], 1.0);
      cE[0][0] = max<double>(cE[0][0], 0.0);
      cE[0][1] = 1.0 - cE[1][0];
      
      indxE[1] = 0;
      
      cE[1][0] = (I1last - I1)/(I1last - I1cur);
      cE[1][0] = min<double>(cE[1][0], 1.0);
      cE[1][0] = max<double>(cE[1][0], 0.0);
      cE[1][1] = 1.0 - cE[0][0];
      
      ok = true;
    }
    
    if (isnan(cE[0][0])) {
      cerr << "NaN in Special Case #1: I1last=" << I1last
	   << "  I1min=" << I1X[indxX][numE[0]-minE-1]
	   << "  I1pen=" << I1X[indxX+1][numE[0]-minE]
	   << "  I1int=" << odd2(EX[indxX+1][0], EX[indxX], I1X[indxX])
	   << "  I1=" << I1
	   << "  dI1top=" << I1-I1last
	   << "  dI1bot=" << I1cur-I1last
	   << "  cX[0]=" << cX[0] << "  CX[1]=" << cX[1]
	   << endl;
    }

    if (isnan(cE[1][0])) {
      cerr << "NaN in Special Case #1: I1cur" << I1cur
	   << "  I1last=" << I1last
	   << "  I1=" << I1
	   << endl;
    }


    if (!ok) {

      noboundary = false;
				// Assume below the grid, 
				// unless we find a match below
      cE[0][0] = cE[1][0] = 1.0;
      cE[0][1] = cE[1][1] = 0.0;
      
      indxE[0] = indxE[1] = 0;
      

      // Special case #2: below 
      //      largest energy endpoint


      for (int i=numE[0]-minE-1; i>=0; i--) {

	I1last = I1cur;
	I1cur = cX[0]*I1X[indxX][i] + cX[1]*I1X[indxX+1][0];
    
	if ( (I1-I1cur)*(I1-I1last) <= 0.0 ) {

	  cE[0][0] = (I1last - I1)/(I1last - I1cur);
	  cE[0][0] = min<double>(cE[0][0], 1.0);
	  cE[0][0] = max<double>(cE[0][0], 0.0);
	  cE[0][1] = 1.0 - cE[1][0];
      
	  ok = true;

	  if (isnan(cE[0][0])) {
	    cerr << "NaN in Special Case #2: I1last=" << I1last
		 << "  I1cur=" << I1cur
		 << "  I1=" << I1
		 << "  cX[0]=" << cX[0]
		 << "  cX[1]=" << cX[1]
		 << "  i=" << i
		 << endl;

	    break;
	  }

	}
      }
      
    }
	
  }
  
  return;
}

bool ResPot::getValues(double I1, double I2,
		       double& O1, double& O2)
{
  O1 = 0.0;
  O2 = 0.0;
  
  int indxX;
  int indxE[2];
  double cX[2];
  double cE[2][2];
  
  bool noboundary;
  getInterp(I1, I2, indxX, indxE, cX, cE, noboundary);
  
  // Interpolate frequencies
  // -----------------------
  
  double fac;
  int num;
  
  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<2; i2++) {
      
      RW *rw = &(orbmat[indxX+i1][indxE[i1]+i2]);
      num = rw->num;
      fac = cX[i1]*cE[i1][i2];
      
      if (ngrid != num) {
	cerr << "ResPot::getValues[1]: Oops! ngrid=" 
	     << ngrid << "  num=" << num 
	     << "  indxX=" << indxX+i1 << "/" << numx
	     << "  indxE=" << indxE[i1]+i2 << "/" << EX[indxX+i1].size()
	     << endl;
      }
      
      O1  += fac * rw->O1;
      O2  += fac * rw->O2;
    }
  }
  
  // return noboundary;
  return true;
}

bool ResPot::getValues(double I1, double I2, CVector& bcoef,
		       double& O1, double& O2,
		       double& Jm, double& dJm,
		       KComplex& Ul, KComplex& dUldE, KComplex& dUldK)
{
  O1 = 0.0;
  O2 = 0.0;
  Jm = 0.0;
  dJm = 0.0;
  Ul = 0.0;
  dUldE = 0.0;
  dUldK = 0.0;
  
  
  int indxX;
  int indxE[2];
  double cX[2];
  double cE[2][2];

  bool wasok = true;

  if (isnan(I1) || isnan(I2)) {
    cerr << "NaN on input values to getInterp\n";
    wasok = false;
  }


  bool noboundary;
  getInterp(I1, I2, indxX, indxE, cX, cE, noboundary);
  
  
  for (int i1=0; i1<2; i1++) {
    if (indxE[i1] >= (int)EX[indxX+i1].size()-1) {
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
      
      RW *rw = &(orbmat[indxX+i1][indxE[i1]+i2]);
      num = rw->num;
      fac = cX[i1]*cE[i1][i2];
      
      if (isnan(fac)) {
	if (wasok) 
	  cerr << "ResPot::getValues[2]: fac=NaN and was OK!!  cX[" << i1 << "]=" << cX[i1]
	       << "  cE[" << i1 << "][" << i2 << "]=" << cE[i1][i2] 
	       << "  indxX=" << indxX+i1 << "/" << numx
	       << "  indxE=" << indxE[i1]+i2 << "/" << EX[i1].size()-1
	       << endl;
	else
	  cerr << "ResPot::getValues[2]: fac=NaN!!  cX[" << i1 << "]=" << cX[i1]
	       << "  cE[" << i1 << "][" << i2 << "]=" << cE[i1][i2] 
	       << "  indxX=" << indxX+i1 << "/" << numx
	       << "  indxE=" << indxE[i1]+i2 << "/" << EX[i1].size()-1
	       << endl;
      }

      if (ngrid != num) {
	assert( EX[indxX+i1].size() == I1X[indxX+i1].size() );
	cerr << "ResPot::getValues[2]: Oops! ngrid=" 
	     << ngrid << "  num=" << num 
	     << "  i1, i2=" << i1 << ", " << i2
	     << "  indxX=" << indxX+i1 << "/" << numx
	     << "  indxE=" << indxE[i1]+i2 << "/" << EX[indxX+i1].size()
	     << endl;
      }
      
      O1  += fac * rw->O1;
      O2  += fac * rw->O2;
      Jm  += fac * rw->Jm;
      dJm += fac * rw->dJm;
      
      for (int n=1; n<=NMAX; n++) {
	Ul +=    fac * rw->W[n-1]   * bcoef[n];
	dUldE += fac * rw->dWE[n-1] * bcoef[n];
	dUldK += fac * rw->dWK[n-1] * bcoef[n];
      }
      
    }
  }
  
  // return noboundary;
  return true;
}


bool ResPot::coord(double* pos, double* vel,
		   double I1, double I2, double beta,
		   double w1, double w2, double w3)
{
  // Linear interpolation coefficients
  // ---------------------------------
  
  int indxX;
  int indxE[2];
  double cX[2];
  double cE[2][2];
  
  bool noboundary;
  getInterp(I1, I2, indxX, indxE, cX, cE, noboundary);
  
  for (int i1=0; i1<2; i1++) {
    if (indxE[i1] >= (int)EX[indxX+i1].size()-1) {
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
  
  double E = 0.0;
  double rmin=Rmax, rmax=Rmin;
  
  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<2; i2++) {
      
      RW *rw = &(orbmat[indxX+i1][indxE[i1]+i2]);
      num = rw->num;
      fac = cX[i1]*cE[i1][i2];
      
      if (ngrid != num) {
	cerr << "ResPot::coord[2]: Oops! ngrid=" 
	     << ngrid << "  num=" << num
	     << "  indxX=" << indxX+i1 << "/" << numx
	     << "  indxE=" << indxE[i1]+i2 << "/" << EX[indxX+i1].size()
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
    vector<double> param(2);
    param[0] = I2;
    param[1] = E;
    
    
    if (w1>2.0*M_PI*0.9 || w1 < 0.1*2.0*M_PI) {
      
      if ( circ.find(get_rc, I2, 0.5*rmin, 2.0*rmax, 1.0e-10, rcirc) ) {
#ifdef DEBUG_VERBOSE
	cout << "  Rcirc bounds error: val, rmin=" << t1 << ", " << rmin
	     << "  val, rmax=" << t2 << ", " << rmax << endl;
#endif
	return false;
      }
      
      if ( tp.find(adj_r, param, Rmin, rcirc, 1.0e-10, r) ) {
#ifdef DEBUG_VERBOSE	
	cout << "  Radj inner bounds error: E=" << E << "  J=" << I2
	     << "  val, r=" << t1 << ", " << Rmin
	     << "  val, rcirc=" << t2 << ", " << rcirc << endl;
#endif
	return false;
      }
      
      rv = 1.0;			// Pericenter
      w1 = 0.0;
      
    } else {
      
      if ( circ.find(get_rc, I2, 0.5*rmin, 2.0*rmax, 1.0e-10, rcirc) ) {
#ifdef DEBUG_VERBOSE
	cout << "  Rcirc outer bounds error: val, rmin=" << t1 << ", " << rmin
	     << "  val, rmax=" << t2 << ", " << rmax << endl;
#endif
	return false;
      }
      
      if ( tp.find(adj_r, param, rcirc, Rmax, 1.0e-10, r) ) {
#ifdef DEBUG_VERBOSE
	cout << "  Radj outer bounds error: E=" << E << "  J=" << I2
	     << "  val, rcirc=" << t1 << ", " << rcirc
	     << "  val, r=" << t2 << ", " << Rmax << endl;
#endif
	return false;
      }
      
      rv = -1.0;		// Apocenter
      w1 = M_PI;
    }
    
    vtot = sqrt(fabs(2.0*(E - halo_model->get_pot(r))));
    vt = I2/r;
    f = 0.0;
    psi = w2;
  }
  
  double vr = 0.0;
  if (vtot > vt) vr = sqrt(vtot*vtot - vt*vt) * rv;
  
  double cost = sinp*sinb;
  double sint = sqrt(fabs(1.0 - cost*cost));
  
  // Check for excursion beyond bounds
  if (fabs(cost)>1.0) {
    cost = copysign(1.0, cost);
    sint = 0.0;
  }
  
  double cos3 = cos(w3);
  double sin3 = sin(w3);
  
  // Check for excursion beyond bounds
  double phi = cosb*cost/(sinb*sint);
  if (fabs(phi)>=1.0) phi = copysign(1.0, phi);
  phi  = asin(phi);
  
  // Phi branch based on Psi
  double tmp  = atan2(sin(psi), cos(psi));
  if (tmp>0.5*M_PI || tmp<-0.5*M_PI) phi = M_PI - phi;
  phi += w3;
  
  double cosf = cos(phi);
  double sinf = sin(phi);
  
  pos[0] = r * sint*cosf;
  pos[1] = r * sint*sinf;
  pos[2] = r * cost;
  
  // Compute velocities
  // ------------------
  
  double xp = cosp*vr - sinp*vt;
  double yp = sinp*vr + cosp*vt;
  
  vel[0] = xp*cos3 - yp*cosb*sin3;
  vel[1] = xp*sin3 + yp*cosb*cos3;
  vel[2] =           yp*sinb;
  
  return noboundary;
}


ofstream* open_debug_file()
{
  ostringstream sout;
  sout << "update." << respot_mpi_id();
  return new ofstream(sout.str().c_str(), ios::app);
}


bool ResPot::Update(double dt, double Phase, double Omega, 
		    double amp, CVector& bcoef,
		    double* posI, double* velI,
		    double* posO, double* velO)
{
  ofstream* out = 0;

  compute_grid();
  
  bool ret;
  
  if (L1==0 && L2==0) {
    for (int k=0; k<3; k++) {
      posO[k] = posI[k];
      velO[k] = velI[k];
    }
    return false;
  }
  
  // Get action angle coords
  //
  double E, K, W1, W2, W3, F, BETA, PSI, I1, I2, O1, O2;
  ret = coord(posI, velI, E, K, I1, I2, O1, O2, W1, W2, W3, F, BETA, PSI);
  if (!ret) return false;
  
  KComplex VB = VeeBeta(L, L2, M, BETA);
  
  double betaM, betaP, beta, BETA0;
  if (BETA-DELTA<0.0) {
    betaM = 0.0;
    betaP = DELTA;
    beta = 0.5*DELTA;
  } else if (BETA+DELTA>M_PI) {
    betaM = M_PI - DELTA;
    betaP = M_PI;
    beta = M_PI - 0.5*DELTA;
  } else {
    betaM = BETA - DELTA;
    betaP = BETA + DELTA;
    beta = BETA;
  }
  BETA0 = BETA;
  
  KComplex DVB = 
    (VeeBeta(L, L2, M, betaP) - VeeBeta(L, L2, M, betaM)) / (betaP - betaM);
  
  // Iterative implicit solution
  //
  double Is0, If0, Is1=0.0, Is2, Is, Is3=0.0;
  double ws0, wf0, ws1=0.0, ws2, ws, ws3=0.0;
  double wg0, wg2, Ig0, Ig2;
  double Jm, dJm, dEIs=0.0, dKIs=0.0;
  KComplex Fw, FI, Ul, Fg, dUldE, dUldK, UldVdIs;
  bool done = false;
  
#ifdef DEBUG_VERBOSE
  double I10 = I1;
  double I20 = I2;
#endif
  
  if (isnan(I1) || isnan(I2)) {
    cerr << "Have a cow!\n";
  }

  // Transformation to slow-fast variables
  //

  ws0 = W1*L1 + W2*L2 + (W3 - Phase)*M;
    
  if (L2) {
    Is2 = Is0 = I2/L2;
    If0 = I1 - Is0*L1;
    wf0 = W1;
  } else {
    Is2 = Is0 = I1/L1;
    If0 = I2 - Is0*L2;
    wf0 = W2;
  }
  wg0 = W3;
  Ig0 = I2*cos(BETA) - Is0*M;
  
  // Angle "drift" for 2nd order calculation
  // 
  if (second_order) {
    ws0 += (O1*L1 + O2*L2 - Omega*M)*0.5*dt;
    if (L2)
      wf0 += O1*0.5*dt;
    else
      wf0 += O2*0.5*dt;
  }
  ws2 = ws0;

  int i;
  for (i=0; i<ITMAX; i++) {
    
    // For convergence test
    //
    ws3 = ws1;
    Is3 = Is1;
    
    // Save previous step
    //
    ws1 = ws2;
    Is1 = Is2;
    
    // For force evaluation
    //
    if (second_order) {
      Is = 0.5*(Is2 + Is0);
      ws = 0.5*(ws2 + ws0);
    } else {
      Is = Is2;
      ws = ws2;
    }
    
    // Canonical transform
    //
    if (L2) {
      I1 = If0 + Is*L1;
      I2 = Is*L2;
    } else {
      I1 = Is*L1;
      I2 = If0 + Is*L2;
    }
    
    getValues(I1, I2, bcoef, O1, O2, Jm, dJm, Ul, dUldE, dUldK);

    // Sanity check
    //
    if (isnan(I1) || isnan(I2)) {
      cerr << "I1 or I2 is NaN: Is0=" 
	   << Is0 << " Is=" << Is << " If0=" 
	   << If0 << " Is2=" << Is2 << " i=" 
	   << i << endl;

      pthread_mutex_lock(&iolock);
      out = open_debug_file();
      *out <<  "I1 or I2 is NaN: Is0=" 
	   << Is0 << " Is=" << Is << " If0=" 
	   << If0 << " Is2=" << Is2 << " i=" 
	   << i << endl;

      out->close();
      pthread_mutex_unlock(&iolock);
      out = 0;
    }


    UldVdIs = Ul * DVB * amp * L2/I2 / tan(beta);
    Ul *= VB * amp;
    dUldE *= VB * amp;
    dUldK *= VB * amp;
    
    dEIs = O1*L1 + O2*L2;
    dKIs = 1.0/Jm*L2;
    
    Fw = (dUldE*dEIs + dUldK*dKIs + UldVdIs)*exp(I*ws);
    
    FI = -I*Ul*exp(I*ws);
    
    // Sanity check
    //
    if (isnan(Fw.real()) || isnan(FI.real())) {
      cerr << "Fw or FI is NaN, dJm=" << dJm 
	   << " Ul="	<< Ul 
	   << " dUldE="	<< dUldE 
	   << " dUldK="	<< dUldK 
	   << " dEIs="	<< dEIs 
	   << " dKIs="	<< dKIs 
	   << " O1="	<< O1 
	   << " O2="	<< O2 
	   << " Omega="	<< Omega
	   << " dt="	<< dt
	   << " ws="	<< ws
	   << " ws0="	<< ws0
	   << " ws2="	<< ws2
	   << " i="	<< i << endl;

      pthread_mutex_lock(&iolock);
      out = open_debug_file();
      *out  << "Fw or FI is NaN, dJm=" << dJm 
	   << " Ul="	<< Ul 
	   << " dUldE="	<< dUldE 
	   << " dUldK="	<< dUldK 
	   << " dEIs="	<< dEIs 
	   << " dKIs="	<< dKIs 
	   << " O1="	<< O1 
	   << " O2="	<< O2 
	   << " Omega="	<< Omega
	   << " dt="	<< dt
	   << " ws="	<< ws
	   << " ws0="	<< ws0
	   << " ws2="	<< ws2
	   << " i="	<< i << endl;
      out->close();
      pthread_mutex_unlock(&iolock);
      out = 0;
    }

    // Update
    
    ws2 = ws0 + dt*Fw.real();
    Is2 = Is0 + dt*FI.real();
    
    if (isnan(ws2)) {
      cerr << "ws2 is NaN, Fw=" << Fw.real()
	   << " ws0=" << ws0
	   << " dt=" << dt
	   << " i="	<< i << endl;

      pthread_mutex_lock(&iolock);
      out = open_debug_file();
      *out  << "ws2 is NaN, Fw=" << Fw.real()
	   << " ws0=" << ws0
	   << " dt=" << dt
	   << " i="	<< i << endl;
      out->close();
      pthread_mutex_unlock(&iolock);
      out = 0;
    }

    
    // Check for convergence
    
    if (fabs(ws1-ws2)<TOLITR && 
	fabs(Is1-Is2)/(fabs(Is2)+1.0e-10)<TOLITR) done = true;
    
    // Limit 2-cycle detection
    
    if (i>3 &&
	fabs(ws3-ws2)<TOLITR && 
	fabs(Is3-Is2)/(fabs(Is2)+1.0e-10)<TOLITR) done = true;
    
    if (done) break;
    
  }
  
  // Update orbital plane angle

  Fg = -Ul * DVB * amp * exp(I*ws)/(sin(beta)*I2);
  wg2 = wg0 + dt*Fg.real();
  Ig2 = Ig0;

  if (!done) {
    // #ifdef DEBUG
    cerr << "Phase, E, K, I1, I2, DI, Dw, Ul, dUldE, dUldK, dEIs, dKIs = " 
	 << Phase
	 << ", " << E
	 << ", " << K
	 << ", " << I1
	 << ", " << I2
	 << ", " << Is2-Is1
	 << ", " << ws2-ws1
	 << ", " << Is2-Is3
	 << ", " << ws2-ws3
	 << ", " << Ul
	 << ", " << dUldE
	 << ", " << dUldK
	 << ", " << dEIs
	 << ", " << dKIs
	 << endl;
    // #endif
    ret = 1;
  }
  
  // Update phase
  //

  double afac = 1.0;
  if (second_order) afac = 0.5;

  ws2 += (O1*L1 + O2*L2 - Omega*M)*afac*dt;
  if (L2)
    wf0 += O1*afac*dt;
  else
    wf0 += O2*afac*dt;


  // Canonical transformation from slow-fast to action-angle
  // 
  W3 = wg2;
  if (L2) {
    I1 = If0 + Is2*L1;
    I2 = Is2*L2;
    W1 = wf0;
    W2 = (ws2 - W1*L1 - (W3 - Phase - Omega*dt)*M)/L2;
  } else {
    I1 = Is2*L1;
    I2 = If0 + Is2*L2;
    W2 = wf0;
    W1 = (ws2 - W2*L2 - (W3 - Phase - Omega*dt)*M)/L1;
  }
    
  double I3 = Is2*M + Ig2;
  double cosb = I3/I2;
  cosb = min<double>(cosb,  1.0);
  cosb = max<double>(cosb, -1.0);
  BETA = acos(cosb);
  
#ifdef DEBUG_VERBOSE
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
  //
  ret = coord(posO, velO, I1, I2, BETA, W1, W2, W3); 
  if (!ret) {
    for (int k=0; k<3; k++) {
      posO[k] = posI[k];
      velO[k] = velI[k];
    }
  }
  
  // Debug
  for (int k=0; k<3; k++) {
    if (isnan(posO[k]) || isinf(posO[k]) || isnan(velO[k]) || isinf(velO[k]))
      ret = 0;
  }
  
  
  return ret;
}


bool ResPot::Force(double dt, double Phase, double Omega, 
		   double amp, CVector& bcoef,
		   double* pos, double* vel, double* acc)
{
  bool ret;
  double pos0[3], vel0[3], pos2[3], vel2[3];
  double E, K, W1, W2, W3, F, BETA, PSI, I1, I2, O1, O2;
  
  // Get action angle coords
  ret = coord(pos, vel, E, K, I1, I2, O1, O2, W1, W2, W3, F, BETA, PSI);
  
  if (!ret) {
    for (int k=0; k<3; k++) acc[k] = 0.0;
    return ret;
  }
  
  // Get phase space update without perturbation
  W1 += O1*dt;
  W2 += O2*dt;
  ret = coord(pos0, vel0, I1, I2, BETA, W1, W2, W3);
  
  if (!ret) {
    for (int k=0; k<3; k++) acc[k] = 0.0;
    return ret;
  }
  
  // Get phase space update with perturbation
  ret = Update(dt, Phase, Omega, amp, bcoef, pos, vel, pos2, vel2);
  
  if (!ret) {
    for (int k=0; k<3; k++) acc[k] = 0.0;
    return ret;
  }
  
  // Effective acceleration
  for (int k=0; k<3; k++) acc[k] = (vel2[k] - vel0[k])/dt;
  
  return ret;
}


void ResPot::check_rw(RW& rw)
{
  for (unsigned i=0; i<rw.r.size(); i++)
    if (isnan(rw.r[i])) cout << "RW error: r nan, i=" << i << endl;
  
  for (unsigned i=0; i<rw.w1.size(); i++)
    if (isnan(rw.w1[i])) cout << "RW error: w1 nan, i=" << i << endl;
  
  for (unsigned i=0; i<rw.f.size(); i++)
    if (isnan(rw.f[i])) cout << "RW error: f nan, i=" << i << endl;
  
  for (unsigned i=0; i<rw.W.size(); i++)
    if (isnan(rw.W[i])) cout << "RW error: W nan, i=" << i << endl;
  
  for (unsigned i=0; i<rw.dWE.size(); i++)
    if (isnan(rw.dWE[i])) cout << "RW error: dWE nan, i=" << i << endl;
  
  for (unsigned i=0; i<rw.dWK.size(); i++)
    if (isnan(rw.dWK[i])) cout << "RW error: dWK nan, i=" << i << endl;
  
  if (isnan(rw.O1))	cout << "RW error: O1 nan" << endl;
  if (isnan(rw.O2))	cout << "RW error: O2 nan" << endl;
  if (isnan(rw.Jm))	cout << "RW error: Jm nan" << endl;
  if (isnan(rw.dJm))	cout << "RW error: dJm nan" << endl;
  if (isnan(rw.E))	cout << "RW error: E nan" << endl;
  if (isnan(rw.K))	cout << "RW error: K nan" << endl;
  if (isnan(rw.I1))	cout << "RW error: I1 nan" << endl;
}

