#include <iostream>
#include <fstream>
#include <iomanip>

#include <values.h>
#include <ResPot.H>
#include <localmpi.h>

// #define TDEBUG
#undef TDEBUG
#ifdef TDEBUG
#include <pthread.h>  
pthread_mutex_t iolock = PTHREAD_MUTEX_INITIALIZER;
const double E0=-1.64716;
const double K0= 0.84168;
bool first = true;
bool start = true;
#endif

double ResPot::ALPHA = 0.33;
double ResPot::DELTA = 0.01;
double TOLITR = 1.0e-8;
int ResPot::NREC = 40;
int ResPot::NUME = 100;
int ResPot::NUMK = 20;
int ResPot::NUMI = 1000;
int ResPot::ITMAX = 100;
KComplex ResPot::I(0.0, 1.0);


int level_surface(Matrix& z, double val, Vector& xr, Vector& yr, int& num);

EKinterp::EKinterp() :
  emin(0.0), emax(0.0), kmin(0.0), kmax(0.0), nume(2), numk(2)
{
  setsize(1, nume, 1, numk);
}

EKinterp::EKinterp(double Emin, double Emax, double Kmin, double Kmax,
		   int numE, int numK) : 
  emin(Emin), emax(Emax), kmin(Kmin), kmax(Kmax), nume(numE), numk(numK)
{
  setsize(1, nume, 1, numk);
}

void EKinterp::reset(double Emin, double Emax, double Kmin, double Kmax,
		     int numE, int numK)
{
  emin = Emin;
  emax = Emax;
  kmin = Kmin;
  kmax = Kmax;
  nume = numE;
  numk = numK;

  dE = (Emax - Emin)/(numE-1);
  dK = (Kmax - Kmin)/(numK-1);

  setsize(1, nume, 1, numk);
}

double EKinterp::value(double E, double K)
{
  int indx = (int)( (E-emin)/dE );
  int indy = (int)( (K-kmin)/dK );

  indx = max<int>(indx, 0);
  indx = min<int>(indx, nume-2);

  indy = max<int>(indy, 0);
  indy = min<int>(indy, numk-2);

  double a[2], b[2];

  a[0] = (emin+(indx+1)*dE - E)/dE;
  b[0] = 1.0 - a[0];
  
  a[1] = (kmin+(indy+1)*dK - K)/dK;
  b[1] = 1.0 - a[1];
  
  return 
    a[0]*a[1]*(*this)[indx+1][indy+1] +
    b[0]*a[1]*(*this)[indx+2][indy+1] +
    a[0]*b[1]*(*this)[indx+1][indy+2] +
    b[0]*b[1]*(*this)[indx+2][indy+2] ;
}


double EKinterp::deriv(double E, double K)
{
  int indx = (int)( (E-emin)/dE );
  int indy = (int)( (K-kmin)/dK );

  indx = max<int>(indx, 0);
  indx = min<int>(indx, nume-2);

  indy = max<int>(indy, 0);
  indy = min<int>(indy, numk-2);

  double a[2], b[2];

  a[0] = -1.0/dE;
  b[0] =  1.0/dE;
  
  a[1] = -1.0/dK;
  b[1] =  1.0/dK;
  
  return 
    a[0]*a[1]*(*this)[indx+1][indy+1] +
    b[0]*a[1]*(*this)[indx+2][indy+1] +
    a[0]*b[1]*(*this)[indx+1][indy+2] +
    b[0]*b[1]*(*this)[indx+2][indy+2] ;
}



void EKinterp::eval(double E, double K, double& value, double& deriv)
{
  int indx = (int)( (E-emin)/dE );
  int indy = (int)( (K-kmin)/dK );

  indx = max<int>(indx, 0);
  indx = min<int>(indx, nume-2);

  indy = max<int>(indy, 0);
  indy = min<int>(indy, numk-2);

  double a[2], b[2];

  a[0] = (emin+(indx+1)*dE - E)/dE;
  b[0] = 1.0 - a[0];
  
  a[1] = (kmin+(indy+1)*dK - K)/dK;
  b[1] = 1.0 - a[1];
  

  value = 
    a[0]*a[1]*(*this)[indx+1][indy+1] +
    b[0]*a[1]*(*this)[indx+2][indy+1] +
    a[0]*b[1]*(*this)[indx+1][indy+2] +
    b[0]*b[1]*(*this)[indx+2][indy+2] ;

  a[0] = -1.0/dE;
  b[0] =  1.0/dE;

  a[1] = -1.0/dK;
  b[1] =  1.0/dK;
  
  deriv = 
    a[0]*a[1]*(*this)[indx+1][indy+1] +
    b[0]*a[1]*(*this)[indx+2][indy+1] +
    a[0]*b[1]*(*this)[indx+1][indy+2] +
    b[0]*b[1]*(*this)[indx+2][indy+2] ;
}


ResPot::ResPot(AxiSymModel *mod, AxiSymBiorth *bio, 
	       int l, int m, int l1, int l2, int nmax)
{
  halo_model = mod;
  halo_ortho = bio;

  grid_computed = false;
  actions_computed = false;

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
  orb->set_biorth(*halo_ortho, L, NMAX);

  compute_grid();
  compute_actions();
}

ResPot::~ResPot()
{
  delete orb;
}

void ResPot::compute_grid() 
{
  if (grid_computed) return;

  double E, K;
  Vector t(1, NMAX);

  dE = (Emax - Emin)/(NUME-1);
  dK = (Kmax - Kmin)/(NUMK-1);

#ifdef TDEBUG
  ofstream tout;
  if (myid==0) tout.open("respot.chk");
#endif

  for (int i=0; i<NUME; i++) {
    
    E = Emin + dE*i;

    EE.push_back(E);

    ovector orbvec;

    for (int k=0; k<NUMK; k++) {
      RW rw;

      K = Kmin + dK*k;

      if (i==0) KK.push_back(K);

      orb->new_orbit(E, K);
      orb->set_biorth(*halo_ortho, L, NMAX);
      orb->pot_trans(L1, L2, t);
      for (int n=0; n<NMAX; n++) rw.W.push_back(t[n+1]);

      struct ANGLE_GRID * grid = orb->get_angle_grid();

#ifdef TDEBUG
      if (myid=0) {
	tout << setw(15) << E 
	     << setw(15) << K
	     << setw(15) << grid->t[1][0]
	     << setw(15) << grid->r[1][0]
	     << setw(15) << grid->w1[1][0]
	     << setw(15) << grid->f[1][0]
	     << endl;
	if (k==NUMK-1) tout << endl;
      }
#endif
	
      for (int j=0; j<grid->num; j++) {
	rw.r.push_back(grid->r[1][j]);
	rw.w1.push_back(grid->w1[1][j]);
	rw.f.push_back(grid->f[1][j]);
      }
      rw.num = grid->num;

      orbvec.push_back(rw);
      
    }

    Jmax.push_back(orb->Jmax());
    orbmat.push_back(orbvec);
  }

  ngrid = orb->get_angle_grid()->num;

#ifdef TDEBUG
  if (myid==0) tout.close();
#endif

  grid_computed = true;
}


int ResPot::coord(double* ps1, double* vel,
		  double& E, double& K, double& J,
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
  if (fabs(r2)>1.0e-8 && fabs(v2)>1.0e-8) rv /= sqrt(r2*v2);

  if (r>halo_model->get_max_radius()) return 0;

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
  
  
  // Linear interpolation coefficients
  // ---------------------------------

  double cE[2], cK[2], cEd[2], cKd[2];

  E = max<double>(E, Emin);
  E = min<double>(E, Emax);

  int indxE = (int)( (E-Emin)/dE );

  indxE = max<int>(indxE, 0);
  indxE = min<int>(indxE, NUME-2);

  cE[0] = (EE[indxE+1] - E)/dE;
  cE[1] = 1.0 - cE[0];

  cEd[0] = -1.0/dE;
  cEd[1] =  1.0/dE;
    
  double Jm  =  cE[0]*Jmax[indxE] +  cE[1]*Jmax[indxE+1];
    
  K = J/Jm;
  K = max<double>(K, Kmin);
  K = min<double>(K, Kmax);

  int indxK = (int)( (K-Kmin)/dK );

  indxK = max<int>(indxK, 0);
  indxK = min<int>(indxK, NUMK-2);
    
  cK[0] = (KK[indxK+1] - K)/dK;
  cK[1] = 1.0 - cK[0];

  cKd[0] = -1.0/dK;
  cKd[1] =  1.0/dK;


  // Compute angles
  // --------------
  
  double fac;
  int num;

  vector<double> tw(ngrid, 0.0);
  vector<double> tf(ngrid, 0.0);
  vector<double> tr(ngrid, 0.0);

  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<2; i2++) {
      RW *rw = &(orbmat[indxE+i1][indxK+i2]);
      num = rw->num;
      fac = cE[i1]*cK[i2];
      
      if (ngrid != num) {
	cerr << "Oops! ngrid=" << ngrid << "  num=" << num << endl;
      }

      for (int k=0; k<ngrid; k++) {
	tw[k] += fac * rw->w1[k];
	tf[k] += fac * rw->f[k];
	tr[k] += fac * rw->r[k];
      }

    }
  }

#ifdef TDEBUG
				// Print the interpolation array 
				// on the first time through . . . 
  if (first && fabs(E-E0)<0.01 && fabs(K-K0)<0.01) {
    pthread_mutex_lock(&iolock);
    cout << "TDEBUG: Lock acquired\n";
    if (first) {		// Being very cautions . . . 
      cout << "TDEBUG: first\n";
      ofstream tout("respot.tst");
      if (tout) cout << "TDEBUG: file is ok . . . writing\n";
      
      tout << "# ngrid=" << ngrid 
	   << "  rw_ngrid=" << orbmat[indxE][indxK].num << endl;
      
      for (int n=0; n<ngrid; n++) {
	tout << setw(15) << tw[n]
	     << setw(15) << tf[n]
	     << setw(15) << tr[n];
	for (int i1=0; i1<2; i1++) {
	  for (int i2=0; i2<2; i2++) {
	    RW *rw = &(orbmat[indxE+i1][indxK+i2]);
	    tout << setw(15) << rw->w1[n]
		 << setw(15) << rw->f[n]
		 << setw(15) << rw->r[n];
	  }
	}
	tout << endl;
      }
      first = false;
    } else {
      cout << "TDEBUG: you beat me to it\n";
    }
    pthread_mutex_unlock(&iolock);
  }
#endif


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
  
  W1  = w1;
  W2  = w2;
  W3  = w3;
  F   = f;
  PSI = psi;

  return 1;
}


//
// Coordindate mapping to unit interval
//

double ResPot::xJ(double J, double Jmin, double Jmax)
{
  if (J>Jmax) return 1.0;
  if (J<Jmin) return 0.0;
  return pow( (J-Jmin)/(Jmax-Jmin), ALPHA );
}

double ResPot::Jx(double x, double Jmin, double Jmax)
{
  if (x>1.0) return Jmax;
  if (x<0.0) return Jmin;
  return Jmin + (Jmax - Jmin)*pow(x, 1.0/ALPHA);
}


double ResPot::dxJ(double J, double Jmin, double Jmax)
{
  if (J>=Jmax) return 1.0/(Jmax-Jmin);
  if (J<=0.0) return MAXFLOAT;
  return pow( (J-Jmin)/(Jmax-Jmin), ALPHA-1.0 )/(Jmax - Jmin);
}


void ResPot::compute_actions()
{

  if (actions_computed) return;

  compute_grid();

  double E, K;
  I1min =  1.0e30;
  I1max = -1.0e30;
  
  I1_EK.reset(Emin, Emax, Kmin, Kmax, NUME, NUMK);
  O1_EK.reset(Emin, Emax, Kmin, Kmax, NUME, NUMK);
  O2_EK.reset(Emin, Emax, Kmin, Kmax, NUME, NUMK);
  
  for (int i=1; i<=NUME; i++) {
    E = Emin + dE*(i-1);

    for (int j=1; j<=NUMK; j++) {
      K = Kmin + dK*(j-1);

      orb->new_orbit(E, K);
      I1_EK[i][j] = orb->get_action(1);
      O1_EK[i][j] = orb->get_freq(1);
      O2_EK[i][j] = orb->get_freq(2);

      I1min = min<double>(I1min, I1_EK[i][j]);
      I1max = max<double>(I1max, I1_EK[i][j]);
    }
  }

  
#ifdef DEBUG
  cout << "I1min=" << I1min << endl;
  cout << "I1max=" << I1max << endl;
#endif

  int num, nseg;
  Vector ER, KR;
  double X, I1;
  double delE = Emax - Emin;
  double delK = Kmax - Kmin;

  // Let   x = [(I1 - I1_min)/(I1_max - I1_min)]^alpha
  // Then  I1 = I1_min + (I1_max - I1_min)] * x^(1/alpha)

  for (int i=0; i<NUMI; i++) {
    X = (0.5+i)/NUMI;
    I1 = Jx(X, I1min, I1max);
    nseg = level_surface(I1_EK, I1, ER, KR, num);

#ifdef DEBUG
    cout << setw(5) << i << setw(15) << I1 
	 << setw(5) << nseg << setw(5) << num << endl;
#endif

    if (nseg==0 || num<2) continue;

    vector<double> ee, jj;
    for (int j=0; j<num; j++) {
      E = Emin + delE*ER[j+1];
      K = Kmin + delK*KR[j+1];
      orb->new_orbit(E, K);
      ee.push_back(E);
      jj.push_back(K*orb->Jmax());
    }

    I1X.push_back(X);
    I1E.push_back(ee);
    I1J.push_back(jj);
  }

  actions_computed = true;

}


void ResPot::getPert(double I1, double I2, CVector& bcoef,
		     double& Jm, double& dJm,
		     KComplex& Ul, KComplex& dUldE, KComplex& dUldK)
{
  Ul = 0.0;
  dUldE = 0.0;
  dUldK = 0.0;

  // Linear interpolation coefficients
  // ---------------------------------

  double cE[2], cK[2], cEd[2], cKd[2];

  I1 = max<double>(I1min, I1);
  I1 = min<double>(I1max, I1);

  double E = getE(I1, I2);
  E = max<double>(E, Emin);
  E = min<double>(E, Emax);

  int indxE = (int)( (E-Emin)/dE );

  indxE = max<int>(indxE, 0);
  indxE = min<int>(indxE, NUME-2);

  cE[0] = (EE[indxE+1] - E)/dE;
  cE[1] = 1.0 - cE[0];

  cEd[0] = -1.0/dE;
  cEd[1] =  1.0/dE;
    
  Jm  =  cE[0]*Jmax[indxE] +  cE[1]*Jmax[indxE+1];
  dJm = cEd[0]*Jmax[indxE] + cEd[1]*Jmax[indxE+1];
    
  double K = I2/Jm;
  K = max<double>(K, Kmin);
  K = min<double>(K, Kmax);

  int indxK = (int)( (K-Kmin)/dK );

  indxK = max<int>(indxK, 0);
  indxK = min<int>(indxK, NUMK-2);
    
  cK[0] = (KK[indxK+1] - K)/dK;
  cK[1] = 1.0 - cK[0];

  cKd[0] = -1.0/dK;
  cKd[1] =  1.0/dK;


  // Compute Ul
  // ----------
  
  double fac;
  int num;

  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<2; i2++) {
      RW *rw = &(orbmat[indxE+i1][indxK+i2]);
      num = rw->num;
      fac = cE[i1]*cK[i2];
      
      if (ngrid != num) {
	cerr << "Oops! ngrid=" << ngrid << "  num=" << num << endl;
      }

      for (int n=1; n<=NMAX; n++) {
	Ul += fac * rw->W[n-1] * bcoef[n];
	dUldE += cEd[i1] * cK [i2] * rw->W[n-1] * bcoef[n];
	dUldK += cE [i1] * cKd[i2] * rw->W[n-1] * bcoef[n];
      }
    }
  }

}


double ResPot::getE(double I1, double I2)
{
  compute_actions();

				// Compute closest I1 in grid
  if (I1<I1min || I1>I1max) {
#ifdef DEBUG
    cout << "I1=" << I1 << " is out of range [" 
	 << I1min << ", " << I1max << "],  I2=" << I2 << endl;
#endif
    I1 = max<double>(I1min*(1.0+0.25/NUMI), I1);
    I1 = min<double>(I1max*(1.0-0.25/NUMI), I1);
  }

  int indx;
  double x = xJ(I1, I1min, I1max);

  indx = Vlocate(x, I1X);
  indx = max<int>(0, indx);
  indx = min<int>(I1X.size()-2, indx);

  double cI[2];

  cI[0] = (I1X[indx+1] - x)/(I1X[indx+1] - I1X[indx]);
  cI[1] = 1.0 - cI[0];

  double E = 
    cI[0]*odd2(I2, I1J[indx  ], I1E[indx  ]) + 
    cI[1]*odd2(I2, I1J[indx+1], I1E[indx+1]) ;

  return E;
}

int ResPot::coord(double* pos, double* vel,
	    double I1, double I2, double beta,
	    double w1, double w2, double w3)
{
  // Linear interpolation coefficients
  // ---------------------------------

  double cE[2], cK[2];

  double E = getE(I1, I2);
  E = max<double>(E, Emin);
  E = min<double>(E, Emax);

  int indxE = (int)( (E-Emin)/dE );

  indxE = max<int>(indxE, 0);
  indxE = min<int>(indxE, NUME-2);

  cE[0] = (EE[indxE+1] - E)/dE;
  cE[1] = 1.0 - cE[0];
    
  double Jm = cE[0]*Jmax[indxE] + cE[1]*Jmax[indxE+1];
    
  double J = I2;
  double K = J/Jm;
  K = max<double>(K, Kmin);
  K = min<double>(K, Kmax);

  int indxK = (int)( (K-Kmin)/dK );

  indxK = max<int>(indxK, 0);
  indxK = min<int>(indxK, NUMK-2);

  cK[0] = (KK[indxK+1] - K)/dK;
  cK[1] = 1.0 - cK[0];

  
  // Compute radius, f(w1)
  // --------------------
  
  double fac;
  int num;

  vector<double> tw(ngrid, 0.0);
  vector<double> tf(ngrid, 0.0);
  vector<double> tr(ngrid, 0.0);

  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<2; i2++) {
      RW *rw = &(orbmat[indxE+i1][indxK+i2]);
      num = rw->num;
      fac = cE[i1]*cK[i2];
      
      if (ngrid != num) {
	cerr << "Oops! ngrid=" << ngrid << "  num=" << num << endl;
      }

      for (int k=0; k<ngrid; k++) {
	tw[k] += fac * rw->w1[k];
	tf[k] += fac * rw->f[k];
	tr[k] += fac * rw->r[k];
      }
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
  double vt   = J/r;
  if (vtot < vt) { 
    vt = vtot;
  }
  double vr   = sqrt(fabs(vtot*vtot - vt*vt)) * rv;

  double cost = sinp*sinb;
  double sint = sqrt(1.0 - cost*cost);

  double cos3 = cos(w3);
  double sin3 = sin(w3);

  double phi  = asin(cosb*cost/(sinb*sint));
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

  return 1;
}


void ResPot::Update(double dt, double Phase, double Omega, 
		    double amp, CVector& bcoef,
		    double* posI, double* velI,
		    double* posO, double* velO)
{
  compute_actions();

  if (L1==0 && L2==0) {
    for (int k=0; k<3; k++) {
      posO[k] = posI[k];
      velO[k] = velI[k];
    }
    return;
  }


  double E, K, W1, W2, W3, F, BETA, PSI, I1, I2;

				// Get action angle coords
  coord(posI, velI, E, K, I2, W1, W2, W3, F, BETA, PSI);
  
  KComplex VB = VeeBeta(L, L2, M, BETA);
  

				// Iterative implicit solution
  double Is0, If0, Is1, Is2, Is;
  double ws0, ws1, ws2, ws;
  double Jm, dJm, dEIs, dKIs;
  KComplex Fw, FI, Ul, dUldE, dUldK;
  bool done = false;

  I1 = I1_EK.value(E, K);

  if (L2) {
    Is2 = Is0 = I2/L2;
    If0 = I1 - Is0*L1;
  } else {
    Is2 = Is0 = I1/L1;
    If0 = I2 - Is0*L2;
  }

  ws2 = ws0 = W1*L1 + W2*L2 + (W3 - Phase)*M;
  
  int i;
  for (i=0; i<ITMAX; i++) {

    ws1 = ws2;
    Is1 = Is2;

    // Force

    Is = 0.5*(Is2 + Is0);
    ws = 0.5*(ws2 + ws0);
    
    if (L2) {
      I1 = If0 + Is*L1;
      I2 = Is*L2;
    } else {
      I1 = Is*L1;
      I2 = If0 + Is*L2;
    }

    getPert(I1, I2, bcoef, Jm, dJm, Ul, dUldE, dUldK);

    Ul *= VB * amp;
    dUldE *= VB * amp;
    dUldK *= VB * amp;

    E = getE(I1, I2);
    E = max<double>(E, Emin);
    E = min<double>(E, Emax);

    K = I2/Jm;
    K = max<double>(K, Kmin);
    K = min<double>(K, Kmax);
      
    dEIs = O1_EK.value(E, K)*L1 + O2_EK.value(E, K)*L2;
    dKIs =  1.0/Jm*L2 - K*dJm*dEIs/Jm;

    Fw = O1_EK.value(E, K)*L1 + O2_EK.value(E, K)*L2 - Omega*M + I*Ul*exp(I*ws);
    FI = -(dUldE*dEIs + dUldK*dKIs)*exp(I*ws);


    // Update

    ws2 = ws0 + dt*Fw.real();
    Is2 = Is0 + dt*FI.real();


    // Check for convergence
    
    if (fabs(ws1-ws2)<TOLITR && 
	fabs(Is1-Is2)/(fabs(Is2)+1.0e-10)<TOLITR) done = true;

    if (done) break;

  }

#ifdef DEBUG
  if (!done) {
    cerr << "Convergence error==>  E, I1, I2, DI, Dw = " 
	 << E
	 << ", " << I1
	 << ", " << I2
	 << ", " << Is2-Is1
	 << ", " << ws2-ws1
	 << endl;
  }
#endif

  // Coordinate transform: Action-Angle to Cartesian

  if (L2) {
    I1 = If0 + 0.5*(Is2+Is0)*L1;
    I2 = 0.5*(Is2+Is0)*L2;
  } else {
    I1 = 0.5*(Is2+Is0)*L1;
    I2 = If0 + 0.5*(Is2+Is0)*L2;
  }

  E = getE(I1, I2);
  E = max<double>(E, Emin);
  E = min<double>(E, Emax);

  K = I2/Jm;
  K = max<double>(K, Kmin);
  K = min<double>(K, Kmax);
      
  if (L2) {
    W1 += O1_EK.value(E, K)*dt;
    W2 = (ws2 - W1*L1 - (W3 - Phase - Omega*dt)*M)/L2;
    I1 = If0 + Is2*L1;
    I2 = Is2*L2;
  } else {
    W2 += O2_EK.value(E, K)*dt;
    W1 = (ws2 - W2*L2 - (W3 - Phase - Omega*dt)*M)/L1;
    I1 = Is2*L1;
    I2 = If0 + Is2*L2;
  }

  coord(posO, velO, I1, I2, BETA, W1, W2, W3);
}


void ResPot::Force(double dt, double Phase, double Omega, 
		   double amp, CVector& bcoef,
		   double* pos, double* vel, double* acc)
{
  double pos2[3], vel2[3];
  
  Update(dt, Phase, Omega, amp, bcoef, pos, vel, pos2, vel2);

  for (int k=0; k<3; k++) acc[k] = (vel2[k] - vel[k])/dt;
}

