#include <values.h>
#include <ResPot.H>
#include <gaussQ.h>

double ResPot::DELTA = 0.01;
int ResPot::NREC = 10;
int ResPot::NUME = 100;
int ResPot::NUMK = 20;
bool ResPot::use_grid = true;
KComplex ResPot::I(0.0, 1.0);

ResPot::ResPot(AxiSymModel *mod, AxiSymBiorth *bio, 
	       int l, int m, int l1, int l2, int nmax)
{
  halo_model = mod;
  halo_ortho = bio;

  grid_computed = false;

  L = l;
  M = m;
  L1 = l1;
  L2 = l2;
  NMAX = nmax;

  Emax = halo_model->get_pot(halo_model->get_max_radius())*(1.0+DELTA);
  Emin = halo_model->get_pot(halo_model->get_min_radius())*(1.0-DELTA);
  Kmin = DELTA;
  Kmax = 1.0-DELTA;

  lq = new LegeQuad(NREC);

  // SphericalOrbit::RMAXF=1.0;
  orb = new SphericalOrbit(halo_model);
  orb->set_biorth(*halo_ortho, L, NMAX);
}

ResPot::~ResPot()
{
  delete lq;
  delete orb;
}

void ResPot::compute_grid() 
{
  if (grid_computed) return;

  double E, K;
  Vector t(1, NMAX);

  dE = (Emax - Emin)/(NUME-1);
  dK = (Kmax - Kmin)/(NUMK-1);

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
  tw = vector<double>(ngrid);
  tf = vector<double>(ngrid);
  tr = vector<double>(ngrid);

  grid_computed = true;
}

double ResPot::Pot(double* pos, double* vel, double phase, CVector &bcoef)
{
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

  if (r>halo_model->get_max_radius()) return 0.0;

  if (use_grid) compute_grid();


  // Compute E, J, beta
  // ------------------

  double E = 0.5*v2 + halo_model->get_pot(r);
  
  double angmom[3];
  
  angmom[0] = pos[1]*vel[2] - pos[2]*vel[1];
  angmom[1] = pos[2]*vel[0] - pos[0]*vel[2];
  angmom[2] = pos[0]*vel[1] - pos[1]*vel[0];

  double J = 0.0, beta = 0.0, K;
  for (int i=0; i<3; i++) J += angmom[i]*angmom[i];
  J = sqrt(J);

  if (J>0.0) beta = acos(angmom[2]/J);
  
  if (use_grid) {

    // Linear interpolation coefficients
    // ---------------------------------

    double cE[2], cK[2];

    E = max<double>(E, Emin);
    E = min<double>(E, Emax);

    int indxE = (int)( (E-Emin)/dE );

    indxE = max<int>(indxE, 0);
    indxE = min<int>(indxE, NUME-2);

    cE[0] = (EE[indxE+1] - E)/dE;
    cE[1] = 1.0 - cE[0];
    
    double Jm = cE[0]*Jmax[indxE] + cE[1]*Jmax[indxE+1];
    
    K = J/Jm;
    K = max<double>(K, Kmin);
    K = min<double>(K, Kmax);

    int indxK = (int)( (K-Kmin)/dK );

    indxK = max<int>(indxK, 0);
    indxK = min<int>(indxK, NUMK-2);

    cK[0] = (KK[indxK+1] - K)/dK;
    cK[1] = 1.0 - cK[0];

    double theta=0.0;
    if (r>0.0) theta = acos(pos[2]/r);
    double phi = atan2(pos[1], pos[0]);

    
    // Compute angles
    // --------------
  
    KComplex Ul=0.0;
    double w1=0.0, w2=0.0, w3, f=0.0, fac, psi;
    int num;

    for (int n=0; n<ngrid; n++) tw[n] = tf[n] = tr[n] = 0.0;

    for (int i1=0; i1<2; i1++) {
      for (int i2=0; i2<2; i2++) {
	RW *rw = &(orbmat[indxE+i1][indxK+i2]);
	num = rw->num;
	fac = cE[i1]*cK[i2];

	for (int n=0; n<ngrid; n++) {
	  tw[n] += fac * rw->w1[n];
	  tf[n] += fac * rw->f[n];
	  tr[n] += fac * rw->r[n];
	}

	for (int n=1; n<=NMAX; n++) {
	  Ul += fac * rw->W[n-1] * bcoef[n];
	}
      }
    }

    w1 = odd2(r, tr, tw);
    f  = odd2(r, tr, tf);

    if (rv < 0.0) {
      w1 = 2.0*M_PI - w1;
      f *= -1.0;
    }

    
    w3 = atan2(angmom[1], angmom[0]) + 0.5*M_PI;


    if (fabs(beta)<1.0e-08) psi = phi - w3;
    else {
      double tmp = cos(theta)/sin(beta);

      if (fabs(tmp)>1.0) {
	if (tmp>1.0) psi =  0.5*M_PI;
	else         psi = -0.5*M_PI;
      } 
      else psi = asin(tmp);

				// Map into [-Pi, Pi]
      phi = atan2(sin(phi-w3), cos(phi-w3));

				// Psi branch
      if (phi>0.5*M_PI || phi<-0.5*M_PI) psi = M_PI - psi;
    }

    w2 = psi + f;

    double arga = w1*L1 + w2*L2 + w3*M - phase*M;

    Ul *= VeeBeta(L, L2, M, beta) * exp(I*arga);

    return Ul.real();

  } else {

    E = max<double>(E, Emin);
    E = min<double>(E, Emax);

    orb->new_orbit(E, 0.5);
    
    K = J/orb->Jmax();
    K = max<double>(K, Kmin);
    K = min<double>(K, Kmax);

    orb->new_orbit(E, K);

    double theta=0.0;
    if (r>0.0) theta = acos(pos[2]/r);
    double phi = atan2(pos[1], pos[0]);

    
    // Compute angles
    // --------------
  
    double w1, w2, w3, f, psi;
    orb->get_angle(2, 0.0);
    struct ANGLE_GRID * grid = orb->get_angle_grid();
	
    w1 = odd2(r, grid->r[1], grid->w1[1]);
    f  = odd2(r, grid->r[1], grid->f[1]);
    
    if (rv < 0.0) {
      w1 = 2.0*M_PI - w1;
      f *= -1.0;
    }

    w3 = atan2(angmom[1], angmom[0]) + 0.5*M_PI;

    if (fabs(beta)<1.0e-10) psi = phi - w3;
    else {
      double tmp = cos(theta)/sin(beta);
      
      if (fabs(tmp)>1.0) {
	if (tmp>1.0) psi =  0.5*M_PI;
	else         psi = -0.5*M_PI;
      } 
      else psi = asin(tmp);
    }

				// Map into [-Pi, Pi]
    phi = atan2(sin(phi - w3), cos(phi - w3));

    if (phi>0.5*M_PI || phi<-0.5*M_PI) psi = M_PI - psi;
  
    w2 = psi + f;

    double arga = w1*L1 + w2*L2 + w3*M - phase*M;
    KComplex Ul=0.0, phase=exp(I*arga);
    Vector trns(1, NMAX);
    orb->pot_trans(L1, L2, trns);

    for (int i=1; i<=NMAX; i++) {
      Ul += phase*bcoef[i]*trns[i];
    }
    Ul *= VeeBeta(L, L2, M, beta);

    return Ul.real();
  }
  
}



void ResPot::ForceSph(double* ps1, double* vel, double phase, CVector &bcoef,
		      double& pot, double &fr, double &ft, double &fp)
{
  // Numerical limit if beta=0
  // -------------------------

  double pos[3], rt=0.0;
  for (int i=0; i<3; i++) {
    pos[i] = ps1[i];
    rt += ps1[i]*ps1[i];
  }

  // Particle at origin?
  if (rt < 1.0e-16) {
    pot = Pot(ps1, vel, phase, bcoef);
    fr = ft = fp = 0.0;
    return;
  }

  if (fabs(ps1[2])<1.0e-4*sqrt(rt)) pos[2] = 1.0e-4*sqrt(rt);

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

  if (r>halo_model->get_max_radius()) {
    pot = fr = ft = fp = 0.0;
    return;
  }

  double cost = cos(pos[2]/r);
  double sint = sqrt(1.0 - cost*cost);
  double phi = atan2(pos[1], pos[0]);
  double theta = 0.0;
  if (r>0.0) theta = acos(pos[2]/r);

  if (use_grid) compute_grid();


  // Compute E, J, beta
  // ------------------

  double E = 0.5*v2 + halo_model->get_pot(r);
  
  double angmom[3];
  
  angmom[0] = pos[1]*vel[2] - pos[2]*vel[1];
  angmom[1] = pos[2]*vel[0] - pos[0]*vel[2];
  angmom[2] = pos[0]*vel[1] - pos[1]*vel[0];

  double J = 0.0, beta = 0.0, K, dwr, dfr, crct;
  double w1=0.0, f=0.0, psi;
  for (int i=0; i<3; i++) J += angmom[i]*angmom[i];
  J = sqrt(J);

  if (J>0.0) beta = acos(angmom[2]/J);
  
  KComplex Ul=0.0, dUldE=0.0, dUldK=0.0;

  if (use_grid) {

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
    
    double Jm = cE[0]*Jmax[indxE] + cE[1]*Jmax[indxE+1];
    
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

    dwr = dfr = crct = 0.0;
    for (int n=0; n<ngrid; n++) tw[n] = tf[n] = tr[n] = 0.0;

    for (int i1=0; i1<2; i1++) {
      for (int i2=0; i2<2; i2++) {
	RW *rw = &(orbmat[indxE+i1][indxK+i2]);
	num = rw->num;
	fac = cE[i1]*cK[i2];

	for (int n=0; n<ngrid; n++) {
	  tw[n] += fac * rw->w1[n];
	  tf[n] += fac * rw->f[n];
	  tr[n] += fac * rw->r[n];
	}

	for (int n=1; n<=NMAX; n++) {
	  Ul += fac * rw->W[n-1] * bcoef[n];
	  dUldE += cEd[i1] * cK [i2] * rw->W[n-1] * bcoef[n];
	  dUldK += cE [i1] * cKd[i2] * rw->W[n-1] * bcoef[n];
	}
      }
    }

    w1 = odd2(r, tr, tw);
    f  = odd2(r, tr, tf);
    dwr = drv2(r, tr, tw);
    dfr = drv2(r, tr, tf);

    if (rv < 0.0) {
      w1 = 2.0*M_PI - w1;
      f *= -1.0;
      dwr *= -1.0;
      dfr *= -1.0;
    }

    
    // DEBUG

    /*

    orb->new_orbit(E, K);
    orb->get_angle(2, 0.0);
    struct ANGLE_GRID * grid = orb->get_angle_grid();
	
    double w1t = odd2(r, grid->r[1], grid->w1[1]);
    double ft  = odd2(r, grid->r[1], grid->f[1]);
    double dwrt = drv2(r, grid->r[1], grid->w1[1]);
    double dfrt = drv2(r, grid->r[1], grid->f[1]);
    
    if (rv < 0.0) {
      w1t = 2.0*M_PI - w1t;
      ft *= -1.0;
      dwrt *= -1.0;
      dfrt *= -1.0;
    }

    const double thresh = 0.5;
    double tmp;
    

    if ( (tmp=fabs((w1 - w1t)/w1t)) > thresh ) 
      {
	cout << setw(15) << tmp << setw(15) << w1 << setw(15) << w1t << endl;
      }
    
    if ( (tmp=fabs((f - ft)/ft)) > thresh ) 
      {
	cout << setw(15) << tmp << setw(15) << f << setw(15) << ft << endl;
      }
    
    if ( (tmp=fabs((dwr - dwrt)/dwrt)) > thresh ) 
      {
	cout << setw(15) << tmp << setw(15) << dwr << setw(15) << dwrt << endl;
      }
    
    if ( (tmp=fabs((dfr - dfrt)/dfrt)) > thresh ) 
      {
	cout << setw(15) << tmp << setw(15) << dfr << setw(15) << dfrt << endl;
      }
    
    */

  } else {

    double Ep, Em, delE = 0.003*(Emax - Emin);
    double Kp, Km, delK = 0.01*(Kmax - Kmin);

    E = max<double>(E, Emin);
    E = min<double>(E, Emax);

    orb->new_orbit(E, 0.5);
    
    K = J/orb->Jmax();
    K = max<double>(K, Kmin);
    K = min<double>(K, Kmax);

				// For E derivative
    if (E-0.5*delE<Emin) {
      Ep = Emin + delE;
      Em = Emin;
    } else if (E+0.5*delE>Emax) {
      Ep = Emax;
      Em = Emax - delE;
    } else {
      Ep = E + 0.5*delE;
      Em = E - 0.5*delE;
    }

				// For K derivative
    if (K-0.5*delK<Kmin) {
      Kp = Kmin + delK;
      Km = Kmin;
    } else if (K+0.5*delK>Kmax) {
      Kp = Kmax;
      Km = Kmax - delK;
    } else {
      Kp = K + 0.5*delK;
      Km = K - 0.5*delK;
    }

    orb->new_orbit(E, K);

    double theta=0.0;
    if (r>0.0) theta = acos(pos[2]/r);

    
    // Compute angles
    // --------------
  
    orb->get_angle(2, 0.0);
    struct ANGLE_GRID * grid = orb->get_angle_grid();
	
    w1 = odd2(r, grid->r[1], grid->w1[1]);
    f  = odd2(r, grid->r[1], grid->f[1]);
    dwr = drv2(r, grid->r[1], grid->w1[1]);
    dfr = drv2(r, grid->r[1], grid->f[1]);
    
    if (rv < 0.0) {
      w1 = 2.0*M_PI - w1;
      f *= -1.0;
      dwr *= -1.0;
      dfr *= -1.0;
    }

    Vector trns(1, NMAX);
    orb->pot_trans(L1, L2, trns);

    for (int i=1; i<=NMAX; i++) {
      Ul += bcoef[i]*trns[i];
    }
    
    // Energy deriv
    // ------------
    orb->new_orbit(Ep, K);
    orb->pot_trans(L1, L2, trns);

    for (int i=1; i<=NMAX; i++) {
      dUldE += bcoef[i]*trns[i]/delE;
    }
    
    orb->new_orbit(Em, K);
    orb->pot_trans(L1, L2, trns);

    for (int i=1; i<=NMAX; i++) {
      dUldE += -bcoef[i]*trns[i]/delE;
    }
    
    // Kappa deriv
    // ------------
    orb->new_orbit(E, Kp);
    orb->pot_trans(L1, L2, trns);

    for (int i=1; i<=NMAX; i++) {
      dUldK += bcoef[i]*trns[i]/delK;
    }
    
    orb->new_orbit(E, Km);
    orb->pot_trans(L1, L2, trns);

    for (int i=1; i<=NMAX; i++) {
      dUldK += -bcoef[i]*trns[i]/delK;
    }
    
  }
  

  KComplex VB = VeeBeta(L, L2, M, beta);
  
  // Gradient in beta
  // ----------------
  double betap, betam, db=0.01;;
  if (beta-0.5*db<0.0) {
    betam=0.0;
    betap = db;
  } else if (beta+0.5*db) {
    betam = M_PI-db;
    betap = M_PI;
  } else {
    betam = beta-0.5*db;
    betap = beta+0.5*db;
  }
  KComplex dVBdb = (VeeBeta(L, L2, M, beta)  - VeeBeta(L, L2, M, beta))/db;

  double cosb = cos(beta);
  double sinb = sin(beta);
  double dbdt = -cosb*cosb*cosb*cost/(sint*sint*sint*sinb);
  
  double dpsidt = -(sint + cost*cosb/sinb*dbdt)/sinb/
    sqrt(fabs(1.0 - cost*cost/(sinb*sinb)));

  double dw3dt = (cosb/(sinb*sint*sint) + 
		  cost/(sint*sinb*sinb)*dbdt)/
    sqrt(fabs(1.0 - cost*cost*cosb*cosb/(sint*sint*sinb*sinb)));

  double dKdt =  -K*cosb*cosb*cost/(sint*sint*sint);


  KComplex tmpC;
  
  // Angle computation
  // -----------------

  double w3 = atan2(angmom[1], angmom[0]) + 0.5*M_PI;

  if (fabs(beta)<1.0e-10) psi = phi - w3;
  else {
    double tmp = cos(theta)/sin(beta);
    
    if (fabs(tmp)>1.0) {
      if (tmp>1.0) psi =  0.5*M_PI;
      else         psi = -0.5*M_PI;
    } 
    else psi = asin(tmp);
  }


  // Map Psi into [-Pi, Pi]
  double tmp = atan2(sin(phi - w3), cos(phi - w3));
  if (tmp>0.5*M_PI || tmp<-0.5*M_PI) {
    psi = M_PI - psi;
    dpsidt *= -1.0;
  }
    
  double w2 = psi + f;
  
  double arga = w1*L1 + w2*L2 + w3*M - phase*M;
  KComplex argC = exp(I*arga);

  // Force computation
  // -----------------

  tmpC = Ul * VB * argC;
  pot = tmpC.real();

  tmpC = dUldE*halo_model->get_pot(r)*VB*argC + Ul*VB*argC*(dwr*L1 + dfr*L2);
  fr = -tmpC.real();

  tmpC = Ul * VB * argC * I * M;
  fp = -tmpC.real();
  
  tmpC = dUldK*dKdt*VB*argC + Ul*dVBdb*argC + Ul*VB*argC*(dpsidt*L2 + dw3dt*M);
  ft = -tmpC.real();
  
}


void ResPot::ForceCart(double* ps1, double* vel, double phase, CVector &bcoef,
		       double& pot, double &fx, double &fy, double &fz)
{
  double pos[3];
  for (int i=0; i<3; i++) pos[i] = ps1[i];
  if (fabs(ps1[2])<1.0e-8) pos[2] = 1.0e-8;

  double fr, fp, ft;

  double X = pos[0];
  double Y = pos[1];
  double Z = pos[2];

  double phi = atan2(Y, X);
  double cosp = cos(phi);
  double sinp = sin(phi);

  double cost = Z/sqrt(X*X+Y*Y+Z*Z);
  double sint = sqrt(1.0 - cost*cost);


  ForceSph(pos, vel, phase, bcoef, pot, fr, fp, ft);
  
  fx = fr*sint*cosp + ft*cost*cosp - fp*sinp;
  fy = fr*sint*sinp + ft*cost*sinp + fp*cosp;
  fz = fr*cost      - ft*sint;
}



int ResPot::coord(double* ps1, double* vel, 
		   double& W1, double& W2, double& W3, double& PSI)
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

  if (use_grid) compute_grid();


  // Compute E, J, beta
  // ------------------

  double E = 0.5*v2 + halo_model->get_pot(r);
  
  double angmom[3];
  
  angmom[0] = pos[1]*vel[2] - pos[2]*vel[1];
  angmom[1] = pos[2]*vel[0] - pos[0]*vel[2];
  angmom[2] = pos[0]*vel[1] - pos[1]*vel[0];

  double J = 0.0, beta = 0.0, K;
  for (int i=0; i<3; i++) J += angmom[i]*angmom[i];
  J = sqrt(J);

  if (J>0.0) beta = acos(angmom[2]/J);
  
  E = max<double>(E, Emin);
  E = min<double>(E, Emax);

  orb->new_orbit(E, 0.5);
    
  K = J/orb->Jmax();
  K = max<double>(K, Kmin);
  K = min<double>(K, Kmax);

  orb->new_orbit(E, K);

  double theta=0.0;
  if (r>0.0) theta = acos(pos[2]/r);
  double phi = atan2(pos[1], pos[0]);
  
    
  // Compute angles
  // --------------
  
  double w1, w2, w3, psi, f;
  orb->get_angle(2, 0.0);
  struct ANGLE_GRID * grid = orb->get_angle_grid();
  
  w1 = odd2(r, grid->r[1], grid->w1[1]);
  f  = odd2(r, grid->r[1], grid->f[1]);
    
  if (rv < 0.0) {
    w1 = 2.0*M_PI - w1;
    f *= -1.0;
  }

  w3 = atan2(angmom[1], angmom[0]) + 0.5*M_PI;

  if (fabs(beta)<1.0e-08) psi = phi - w3;
  else {
    double tmp = cos(theta)/sin(beta);

    if (fabs(tmp)>1.0) {
      if (tmp>1.0) psi =  0.5*M_PI;
      else         psi = -0.5*M_PI;
    } 
    else psi = asin(tmp);

				// Map into [-Pi, Pi]
    phi = atan2(sin(phi-w3), cos(phi-w3));

    if (phi>0.5*M_PI || phi<-0.5*M_PI) psi = M_PI - psi;
  
  }

  w2 = psi + f;
  
  W1 = w1;
  W2 = w2;
  W3 = w3;
  PSI = psi;

  return 1;
}
