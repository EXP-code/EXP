#include <iostream>
#include <iomanip>

//#include "shocktube.H"
#include "Initialize.H"
#include "globalInit.H"

double boltz = 1.380648e-16; //in cgs
double mp = 1.672622e-24; //in grams
double amu = 1.660011e-24; //amu


double atomic_masses[2] = {1.00794, 4.002602};


/*void Initialize(std::vector<Particle>& p, double mass,
                double T1, double D1, double T2, double D2,
                vector<double> &L, int nd)
{
  unsigned npart = p.size();
  double   ratio = D1/(D2+D1), v1 = sqrt(T1), v2 = sqrt(T2);
  unsigned first = (unsigned)floor(ratio*npart);
  double rho1    = mass*first/(0.5*L[0]*L[1]*L[2]*npart);
  double rho2    = mass*(npart-first)/(0.5*L[0]*L[1]*L[2]*npart);

  for (unsigned i=0; i<first; i++) {
    p[i].mass = mass/npart;
    for (unsigned k=0; k<3; k++) {
      if (k==0) p[i].pos[k] = 0.5*L[k]*(*Unit)();
      else      p[i].pos[k] = L[k]*(*Unit)();
      p[i].vel[k] = v1*(*Norm)();
    }
    p[i].dattrib.push_back(T1);
    p[i].dattrib.push_back(rho1);
    for (int n=0; n<nd-2; n++) p[i].dattrib.push_back(0.0);
  }

  for (unsigned i=first; i<npart; i++) {
    p[i].mass = mass/npart;
    for (unsigned k=0; k<3; k++) {
      if (k==0) p[i].pos[k] = 0.5*L[k] + 0.5*L[k]*(*Unit)();
      else      p[i].pos[k] = L[k]*(*Unit)();
      p[i].vel[k] = v2*(*Norm)();
    }
    p[i].dattrib.push_back(T2);
    p[i].dattrib.push_back(rho2);
    for (int n=0; n<nd-2; n++) p[i].dattrib.push_back(0.0);
  }
}*/

void InitializeUniform(std::vector<Particle>& p,
                       double mass, double T, vector<double> &L, int nd)
  {
  unsigned npart = p.size();
  //double v0 = sqrt(T);
  double rho = mass/(L[0]*L[1]*L[2]);

  cout << "Temperature: " << T << " K " << endl;
  cout << "Num: " << npart << endl;
  cout << "Length unit: " << Lunit << " cm" << endl;
  cout << "Time unit: " << Tunit << " s" << endl;
  cout << "Vel unit: " << Vunit << " cm/s" << endl;
  cout << "Mass unit: " << Munit << " g" << endl;
  double varH = sqrt((boltz*T)/(atomic_masses[0]*amu));
  double varHe = sqrt((boltz*T)/(atomic_masses[1]*amu));
  for (unsigned i=0; i<npart; i++) {
    p[i].mass = mass/npart;
    //double v0 = sqrt((3*boltz*T)/(atomic_masses[int(p[i].Z)-1]*amu));
    //v0 /= Vunit;
    //cout << v0 << "\t" << endl;

    for (unsigned k=0; k<3; k++) {
      p[i].pos[k] = L[k]*(*Unit)();
      if (p[i].Z == 1) {
     	 //p[i].vel[k] = v0*(*Norm)();
         p[i].vel[k] = varH*(*Norm)();
       }
      if (p[i].Z == 2) p[i].vel[k] = varHe*(*Norm)();
      //p[i].vel[k] = sqrt((boltz*T)/(p[i].mass*Munit))*(*Norm)();
      //cout << fabs(p[i].vel[k]) << "\t";
      p[i].vel[k] /= Vunit;
    }
    p[i].dattrib.push_back(T);
    p[i].dattrib.push_back(rho);
    for (int n=0; n<nd-2; n++) p[i].dattrib.push_back(0.0);
  }
}


/*void InitializeShearingSheet(std::vector<Particle>& p,
                             double mass, double T, vector<double> &L, 
                             double R0, double S0, int nd)
  {
  unsigned npart = p.size();
  double v0 = sqrt(T);
  double rho = mass/(L[0]*L[1]*L[2]);

  double omega = sqrt(2.0)*S0/R0;
  double kappa = 2.0*S0/R0;
  double A     =  S0/(sqrt(2.0)*R0);
  double B     = -S0/(sqrt(2.0)*R0);
  double OmR0  = omega*R0;

  double x, y, u, yd, v, I1, I2, x0, w1;

  for (unsigned i=0; i<npart; i++) {
    p[i].mass = mass/npart;

    for (unsigned k=0; k<3; k++) {
      p[i].pos[k] = L[k]*(*Unit)(); // Uniform
      p[i].vel[k] = v0*(*Norm)();   // From temperature distribution
    }

    x  = p[i].pos[0] - 0.5*L[0];    // Center shear on box center
    y  = p[i].pos[1];
    u  = p[i].vel[0];               // u = xdot
    v  = p[i].vel[1];               // v
    yd = v - 2.0*A*x;               // ydot

    // Compute actions
    //
    I1 = 0.5/kappa*(u*u + 0.25*kappa*kappa/(B*B)*v*v);
    I2 = OmR0 + 2.0*omega*x + yd;

    // Compute w1 angle (and correct branch)
    //
    x0 = (I2 - OmR0)/(-2.0*B);
    w1 = (x - x0)*sqrt(kappa/(2.0*I1));
    if (w1<-1.0)     w1 = -0.5*M_PI;
    else if (w1>1.0) w1 = 0.5*M_PI;
    else             w1 = asin(w1);
    if (u<0.0)       w1 = M_PI - w1;

                                // Sanity check
    double xt = x0 + sqrt(2.0*I1/kappa)*sin(w1);
    if (fabs(xt - x)>1.0e-10) {
      cout << "x=" << x << " xt=" << xt << endl;
    }
                                // xdot
    p[i].vel[0] = sqrt(2.0*kappa*I1)*cos(w1);
                                // ydot
    p[i].vel[1] = 2.0*B/kappa*sqrt(2.0*kappa*I1)*sin(w1) - 2.0*A*x;

    p[i].dattrib.push_back(T);
    p[i].dattrib.push_back(rho);
    for (int n=0; n<nd-2; n++) p[i].dattrib.push_back(0.0);
  }
}


void InitializeUniformWave(std::vector<Particle>& p,
                           double mass, double T, vector<double> &L,
                           double Period, double Amp, int nd)
{
  unsigned npart = p.size();
  double v0 = sqrt(T);

  double delta, deriv, R, x;
  double K = 2.0*M_PI/Period;
  const int MAXIT=40;

  for (unsigned i=0; i<npart; i++) {
    // Mass
    p[i].mass = mass/npart;

    // X position
    x = R = L[0]*(*Unit)();
    for (int j=0; j<MAXIT; j++) {
      delta = x + Amp*sin(K*x)/K - R;
      deriv = 1 + Amp*cos(K*x);
      x += -delta/deriv;
      if (fabs(delta)<1.0e-10) break;
    }
    p[i].pos[0] = x;

    // Y and Z positions
    for (unsigned k=1; k<3; k++) p[i].pos[k] = L[k]*(*Unit)();

    // Velocities
    for (unsigned k=0; k<3; k++) p[i].vel[k] = v0*(*Norm)();

    // Perturbation
    p[i].vel[0] += Amp*v0*cos(K*R);

    // Double attributes
    for (int n=0; n<nd; n++)  p[i].dattrib.push_back(0.0);
  }

}


/*void InitializeUniformShear(std::vector<Particle>& p, double mass,
                            double T1, double D1, double T2, double D2,
                            vector<double> &L, double mach,
                            double Phi, double Period, double Amp,
                            int ni, int nd)
{
  const double gamma = 5/3;     // Adiabatic constant for ideal gas
  unsigned npart = p.size();
  double s1 = sqrt(T1);         // Velocity dispersions
  double s2 = sqrt(T2);
  double c1 = gamma*s1;         // Sound speeds
  double c2 = gamma*s2;

  double vshear = mach*c1;      // Shear speed for Fluid 1
  double K = 2.0*M_PI/Period;   // Perturbation frequency
  double omega1 = Phi*K*c1 - K*vshear;
  double omega2 = Phi*K*c2;
  double vpert1 = Amp*vshear;   // Velocity perturbation amplitudes
  double vpert2 = vpert1*omega1/omega2*D1/D2;

                                // Number fractions for each fluid
  double ratio = D1/(D2+D1);
  unsigned first = (unsigned)floor(ratio*npart);

  double volume = L[0]*L[1]*L[2];
  double epsilon1 = vpert1*omega1/(K*c1*c1);
  double epsilon2 = vpert2*omega2/(K*c2*c2);
  double delta, deriv, R, x;
  const int MAXIT=40;

  if (myid==0) {
    cout << left << endl 
         << setfill('-') << setw(60) << '-' << endl << setfill(' ')
         << "KH initial conditions" << endl
         << setfill('-') << setw(60) << '-' << endl << setfill(' ')
         << setw(10) << ' ' << setw(20) << "Shear" << ": " << vshear << endl
         << setw(10) << ' ' << setw(20) << "K" << ": " << K << endl
         << setw(10) << ' ' << setw(20) << "C_1" << ": " << c1 << endl
         << setw(10) << ' ' << setw(20) << "C_2" << ": " << c2 << endl
         << setw(10) << ' ' << setw(20) << "Omega_1" << ": " << omega1 << endl
         << setw(10) << ' ' << setw(20) << "Omega_2" << ": " << omega2 << endl
         << setw(10) << ' ' << setw(20) << "Vpert_1" << ": " << vpert1 << endl
         << setw(10) << ' ' << setw(20) << "Vpert_2" << ": " << vpert2 << endl
         << setw(10) << ' ' << setw(20) << "Eps_1" << ": " << epsilon1 << endl
         << setw(10) << ' ' << setw(20) << "Eps_2" << ": " << epsilon2 << endl
         << setfill('-') << setw(60) << '-' << endl << setfill(' ') << endl;
  }

  for (unsigned i=0; i<first; i++) {
    // Mass
    p[i].mass = mass/npart;

    // X position
    x = R = L[0]*(*Unit)();
    for (int j=0; j<MAXIT; j++) {
      delta = x + epsilon1*sin(K*x)/K - R;
      deriv = 1 + epsilon1*cos(K*x);
      x += -delta/deriv;
      if (fabs(delta)<1.0e-10) break;
    }
    if (x<0.0)   x += (1.0+floor(-x/L[0]))*L[0];
    if (x>=L[0]) x -= floor(x/L[0])*L[0];
    p[i].pos[0] = x;

    // Y position
    p[i].pos[1] = 0.5*L[1]*(*Unit)();

    // Z position
    p[i].pos[2] = L[2]*(*Unit)();

    // Velocities
    for (unsigned k=0; k<3; k++) p[i].vel[k] = s1*(*Norm)();
    p[i].vel[0] += vshear;

    // Perturbation 
    p[i].vel[0] += vpert1*cos(K*p[i].pos[0]);

    // Attributes
    for (int k=0; k<ni; k++) p[i].iattrib.push_back(0);
    for (int k=0; k<nd; k++) p[i].dattrib.push_back(0.0);
  }

  for (unsigned i=first; i<npart; i++) {
    // Mass
    p[i].mass = mass/npart;

    // X position
    x = R = L[0]*(*Unit)();
    for (int j=0; j<MAXIT; j++) {
      delta = x + epsilon2*sin(K*x)/K - R;
      deriv = 1 + epsilon2*cos(K*x);
      x += -delta/deriv;
      if (fabs(delta)<1.0e-10) break;
    }
    if (x<0.0)   x += (1.0+floor(-x/L[0]))*L[0];
    if (x>=L[0]) x -= floor(x/L[0])*L[0];
    p[i].pos[0] = x;

    // Y position
    p[i].pos[1] = 0.5*L[1] + 0.5*L[1]*(*Unit)();

    // Z position
    p[i].pos[2] = L[2]*(*Unit)();

    // Velocities
    for (unsigned k=0; k<3; k++) p[i].vel[k] = s2*(*Norm)();

    // Perturbation
    p[i].vel[0] += vpert2*cos(K*p[i].pos[0]);
  
    // Attributes
    for (int k=0; k<ni; k++) p[i].iattrib.push_back(0);
    for (int k=0; k<nd; k++) p[i].dattrib.push_back(0.0);
  }

}*/

