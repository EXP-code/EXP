#include <values.h>
#include <ResPot.H>
#include <gaussQ.h>

double ResPot::DELTA = 0.01;
int ResPot::NUME = 400;
int ResPot::NUMK = 100;

KComplex ResPot::I(0.0, 1.0);

ResPot::ResPot(AxiSymModel *mod, AxiSymBiorth *bio, 
	       int l, int m, int l1, int l2, int nmax)
{
  halo_model = mod;
  halo_ortho = bio;

  L = l;
  M = m;
  L1 = l1;
  L2 = l2;
  NMAX = nmax;

  Emax = halo_model->get_pot(halo_model->get_max_radius())*(1.0+DELTA);
  Emin = halo_model->get_pot(halo_model->get_min_radius())*(1.0-DELTA);
  Kmin = DELTA;
  Kmax = 1.0-DELTA;

  // SphericalOrbit::RMAXF=1.0;
  orb = new SphericalOrbit(halo_model);
  orb->set_biorth(*halo_ortho, L, NMAX);

  compute_grid();
}

ResPot::~ResPot()
{
  delete orb;
}

void ResPot::compute_grid() 
{
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

}


double ResPot::Pot(double* pos, double* vel, double Omega, double time,
		   CVector &bcoef)
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
  
  KComplex ans=0.0;
  double w1=0.0, w2=0.0, w3, f=0.0, fac, psi;
  int num;
  
  for (int i1=0; i1<2; i1++) {
    for (int i2=0; i2<2; i2++) {
      RW *rw = &(orbmat[indxE+i1][indxK+i2]);
      num = rw->num;
      fac = cE[i1]*cK[i2];
      
      if (r<rw->r[0]) {
	w1 += fac*rw->w1[0];
	f  += fac*rw->f[0];
      }
      else if (r>rw->r[num-1]) {
	w1 += fac*rw->w1[num-1];
	f  += fac*rw->f[num-1];
      }
      else {
	w1 += fac*odd2(r, rw->r, rw->w1);
	f  += fac*odd2(r, rw->r, rw->f);
      }
      
      for (int n=1; n<=NMAX; n++)
	ans += fac * rw->W[n-1] * bcoef[n];
    }
  }

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

  double arga = w1*L1 + w2*L2 + w3*M - Omega*time*M;
  
  ans *= VeeBeta(L, L2, M, beta) * exp(I*arga);
  
  return ans.real();

}
