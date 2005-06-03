#include <math.h>
#include <pthread.h>

#include <ACG.h>
#include <Uniform.h>
#include <Normal.h>

#include <TwoBodyDiffuse.H>

static pthread_mutex_t randlock = PTHREAD_MUTEX_INITIALIZER;

TwoBodyDiffuse::TwoBodyDiffuse(
			       double PMASS,
			       bool LOGR,
			       double LogL,
			       double Clip,
			       int Numr,
			       int Numv,
			       int Nume,
			       string ModFile,
			       int Diverge,
			       double Diverge_rfac)
{
  pmass = PMASS;		// Particle mass for relaxation
  logr = LOGR;			// Use logarithmic radial intervals
  logL = LogL;			// LogLambda
  clip = Clip;			// Number of sigma for clipping

  numr = 200;			// Number of radial points in model
  numv = 64;			// Number of velocity points in 
				// coef grid
  nume = 128;			// Number of energy integration knots

  modfile = "SLGridSph.model";	// model file rather than histogram
  diverge = Diverge;		// Use power divergence for model
  diverge_rfac = Diverge_rfac;	// with this exponent

  // Set variables

  Gamma = 4.0*M_PI*pmass*logL;
  
  gen = new ACG(seed);
  urand = new Uniform(0.0, 1.0, gen);
  nrand = new Normal(0.0, 1.0, gen);

  compute_diffuse();
}

TwoBodyDiffuse::~TwoBodyDiffuse()
{
  delete model;
  delete jq;
  delete gen;
  delete urand;
  delete nrand;
}

void TwoBodyDiffuse::get_diffusion(double dtime, double* pos, double* vel,
				   double *dvel) 
{
  double rr, vv, vdotr;
  double dvpara1, dvpara2, dvperp2;
  vector<double> er(3), ev(3), e1(3), e2(3);


				// Compute radius and velocity scalars
  rr = vv = vdotr = 0.0;
  for (int k=0; k<3; k++) {
    er[k] = pos[k];
    ev[k] = vel[k];
    rr += er[k]*er[k];
    vv += ev[k]*ev[k];
    vdotr += er[k]*ev[k];
  }

    rr = sqrt(rr+1.0e-12);
    vv = sqrt(vv+1.0e-12);

				// Compute unit vectors
    for (int k=0; k<3; k++) {
      er[k] /= rr;
      ev[k] /= vv;
    }

				// Subtract velocity direction
				// radial direction
    double norm = 0.0;
    for (int k=0; k<3; k++) {
      e1[k] = er[k] - vdotr*ev[k]/(vv*rr);
      norm += e1[k]*e1[k];
    }
    for (int k=0; k<3; k++) e1[k] /= sqrt(norm+1.0e-12);


				// Unit vec perp to ${\hat e1}$ and ${\hat v}$
    e2[0] = ev[1]*e1[2] - ev[2]*e1[1];
    e2[1] = ev[2]*e1[0] - ev[0]*e1[2];
    e2[2] = ev[0]*e1[1] - ev[1]*e1[0];

				// Get diffusion coefficients
    get_coefs(rr, vv, dvpara1, dvpara2, dvperp2);
    

				// Normal random variates
    pthread_mutex_lock(&randlock);
    double r1 = (*nrand)();
    double r2 = (*nrand)();
    double r3 = (*urand)();
    pthread_mutex_unlock(&randlock);

				// Clip

    if (fabs(r1)>clip) r1 = copysign(clip, r1);
    if (fabs(r2)>clip) r2 = copysign(clip, r2);

				// Velocity changes
    double deltaVpara = dvpara1*dtime + r1*sqrt(dvpara2*dtime);
    double deltaVperp = r2*sqrt(dvperp2*dtime);

    double angle = 2.0*M_PI*r3;	// Angle in tangent velocity plane
    double cosA = cos(angle);
    double sinA = sin(angle);
				// Update particle velocity
    for (int k=0; k<3; k++)
      dvel[k] = 
	deltaVpara*ev[k]      +	// Parallel direction
	deltaVperp*e1[k]*cosA +	// Perpendicular direction #1
	deltaVperp*e2[k]*sinA ;	// Perpendicular direction #2

}



void TwoBodyDiffuse::compute_diffuse()
{
  // ===================================================================
  // Create grids (if necessary)
  // =================================================================== 
  
  if (R.size() == 0) {

    R = vector<double>(numr);
    Vmax = vector<double>(numr);
    dVpara1 = vector<dvector>(numr);
    dVpara2 = vector<dvector>(numr);
    dVperp2 = vector<dvector>(numr);

    for (int i=0; i<numr; i++) {
      dVpara1[i] = vector<double>(numv);
      dVpara2[i] = vector<double>(numv);
      dVperp2[i] = vector<double>(numv);
    }

    jq = new LegeQuad(nume);
  }

  model = new SphericalModelTable(modfile, diverge, diverge_rfac);
  model->setup_df(800);

  Rmin = model->get_min_radius();
  Rmax = model->get_max_radius();
  Emin = model->get_pot(Rmin);
  Emax = model->get_pot(Rmax);

  if (logr) {
    delR = (log(Rmax) - log(Rmin))/numr;
    for (int i=0; i<numr; i++) R[i] = Rmin*exp(delR*i);
  }
  else {
    delR = (Rmax - Rmin)/numr;
    for (int i=0; i<numr; i++) R[i] = Rmin + delR*i;
  }


  // ===================================================================
  // Compute diffusion coefficients
  // =================================================================== 
  
  double F0, F1, F2;
  double E, EE, phi, dV, V, df, f1, vrat;

				// Expressions from Cohn 1979
  delV = 1.0/(numv-1);

  for (int i=0; i<numr; i++) {
    
    phi = model->get_pot(R[i]);
    Vmax[i] = sqrt(2.0*(Emax - phi));
    dV = Vmax[i]/numv;

    for (int j=0; j<numv; j++) {

      V = dV*(j+0.5);
      EE = 0.5*V*V + phi;

      F0 = F1 = F2 = 0.0;

      for (int k=1; k<=jq->get_n(); k++) {

	F0 += 2.0*model->distf(EE + (Emax - EE)*jq->knot(k), 0.5) * 
	  jq->weight(k);
      
	E = phi + (EE - phi)*jq->knot(k);
	df = 2.0*model->distf(E, 0.5);
	vrat = (E - phi)/(EE - phi);
	f1 = df*sqrt(vrat);
      
	F1 += f1 * jq->weight(k);
	F2 += f1*vrat * jq->weight(k);
      }
      
      F0 *= 4.0*M_PI*Gamma*(Emax - EE);
      F1 *= 4.0*M_PI*Gamma*(EE - phi);
      F2 *= 4.0*M_PI*Gamma*(EE - phi);
	
      dVpara1[i][j] = -2.0*F1/V;
      dVpara2[i][j] = 2.0*(F0 + F2)/3.0;
      dVperp2[i][j] = 2.0*(2.0*F0 + 3.0*F1 - F2)/3.0;
    
    } // End V loop


  } // End R loop

}


void TwoBodyDiffuse::get_coefs(double r, double v,
			       double& dvpara1, double& dvpara2, double& dvperp2)
{
  
  // ===================================================================
  // Linear interpolation coefficients
  // =================================================================== 

  int indxR, indxV;
  double a0, a1, b0, b1;
  
				// Enforce upper and lower bounds
#ifdef DEBUG
  if (r>1.1*Rmax)
    cout << " Error: R=" << r
	 << " V=" << v
	 << " Rmin=" << Rmin << " Rmax=" << Rmax
	 << ": radius out of bounds" << endl;
#endif      
  r = max<double>(r, Rmin);
  r = min<double>(r, Rmax);

  if (logr)
    indxR = (int)( (log(r) - log(Rmin))/delR );
  else
    indxR = (int)( (r - Rmin)/delR );

  indxR = min<int>(indxR, numr-2);

  if (logr) {
    a0 = (log(R[indxR+1]) - log(r))/delR;
    a1 = (log(r) - log(R[indxR  ]))/delR;
  } else {
    a0 = (R[indxR+1] - r)/delR;
    a1 = (r - R[indxR  ])/delR;
  }

				// Get Vmax
  double vmax = a0*Vmax[indxR] + a1*Vmax[indxR+1];

				// Enforce upper and lower bounds
  double vrel = v/vmax;
#ifdef DEBUG
  if (vrel>1.4)
    cout << "Error: R=" << r
	 << " V=" << v << " Vmax=" << vmax
	 << ": velocity out of bounds" << endl;
#endif      
  vrel = max<double>(vrel, 0.0);
  vrel = min<double>(vrel, 1.0);

  indxV = (int)( vrel/delV-0.5 );
  indxV = min<int>(indxV, numv-2);

  b0 = (delV*(1.5+indxV) - vrel)/delV;
  b1 = (vrel - delV*(0.5+indxV))/delV;


				// Do the interpolation on all quatities
  dvpara1 = 
    b0*(a0*dVpara1[indxR][indxV]   + a1*dVpara1[indxR+1][indxV]  ) +
    b1*(a0*dVpara1[indxR][indxV+1] + a1*dVpara1[indxR+1][indxV+1]) ;

  dvpara2 = 
    b0*(a0*dVpara2[indxR][indxV]   + a1*dVpara2[indxR+1][indxV]  ) +
    b1*(a0*dVpara2[indxR][indxV+1] + a1*dVpara2[indxR+1][indxV+1]) ;

  dvperp2 = 
    b0*(a0*dVperp2[indxR][indxV]   + a1*dVperp2[indxR+1][indxV]  ) +
    b1*(a0*dVperp2[indxR][indxV+1] + a1*dVperp2[indxR+1][indxV+1]) ;

				// Interpolation to unphysical values?
  if (dvpara2 < 0.0 || dvperp2 < 0.0) 
    {
#ifdef DEBUG
      cout << "Error: R=" << r
	   << " V=" << v << " Vmax=" << vmax
	   << " dvpara1, dvpara2, dvperp2="
	   << dvpara1 << ", " << dvpara2 << ", " << dvperp2 << endl;
#endif      
      dvpara1 = dvpara2 = dvperp2 = 0.0;
    }

}

void TwoBodyDiffuse::dump_grid(ostream *out)
{
  out->setf(ios::left);

  *out << setw(15) << "# Radius" 
       << setw(15) << "Velocity"
       << setw(15) << "D(para)"
       << setw(15) << "D(para^2)"
       << setw(15) << "D(perp^2)"
       << endl;

  char c = out->fill('-');
  *out << setw(15) << "#-1" 
       << setw(15) << "|-2"
       << setw(15) << "|-3"
       << setw(15) << "|-4"
       << setw(15) << "|-5"
       << endl;
  out->fill(c);

  double dV;

  for (int i=0; i<numr; i++) {
    dV = Vmax[i]/numv;
    for (int j=0; j<numv; j++) {
      *out << setw(15) << R[i] 
	   << setw(15) << dV*(j+0.5)
	   << setw(15) << dVpara1[i][j]
	   << setw(15) << dVpara2[i][j]
	   << setw(15) << dVperp2[i][j]
	   << endl;
    }
    *out << endl;
  }
  
}
