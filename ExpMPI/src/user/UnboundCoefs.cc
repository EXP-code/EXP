#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

#include <Timer.h>
#include <fft.h>
#include <interp.h>

#include <massmodel.h>
#include <biorth.h>
#include <UnboundCoefs.H>

UnboundCoefs::UnboundCoefs(double Energy, double rperi, 
			   double redge, double delR,
			   AxiSymModel *model, string OUTFILE) :
  E(Energy), Rperi(rperi), Redge(redge), deltaR(delR)
{
  timer.Microseconds();

  // =================================
  // Compute orbit
  // =================================

  double VTperi = sqrt(2.0*(E - model->get_pot(Rperi)));
  double J = Rperi*VTperi;

  R.push_back(Rperi);
  T.push_back(0.0);
  PHI.push_back(0.0);

  //
  // Trapezoidal increments
  //

  double rnext, rlast = Rperi;
  double tnext, tlast = 0.0;
  double phinext, philast = 0.0;

  // First step

  double denom = sqrt(2.0*(VTperi*VTperi/Rperi - model->get_dpot(Rperi)));

  rnext = rlast + deltaR;
  tnext = tlast + 2.0*sqrt(rnext - Rperi)/denom;
  phinext = philast + 2.0*sqrt(rnext - Rperi)/denom * J/(Rperi*Rperi);
  
  R.push_back(rnext);
  T.push_back(tnext);
  PHI.push_back(phinext);

  rlast = rnext;
  tlast = tnext;
  philast = phinext;

  while (R.back() < Redge) {
    rnext = rlast + deltaR;
    tnext = tlast + 0.5*(rnext - rlast)*
      (
       1.0/sqrt(2.0*(E - model->get_pot(rlast)) - J*J/(rlast*rlast)) +
       1.0/sqrt(2.0*(E - model->get_pot(rnext)) - J*J/(rnext*rnext)) 
       );
    phinext = philast + 0.5*(rnext - rlast)*
      (
       J/(rlast*rlast) /
       sqrt(2.0*(E - model->get_pot(rlast)) - J*J/(rlast*rlast)) +
       J/(rnext*rnext) /
       sqrt(2.0*(E - model->get_pot(rnext)) - J*J/(rnext*rnext)) 
       );

    rlast = rnext;
    tlast = tnext;
    philast = phinext;

    R.push_back(rnext);
    T.push_back(tnext);
    PHI.push_back(phinext);
  }
  
  tmax = T.back();

  if (OUTFILE.length()) {
    string outfile = OUTFILE + ".orbit";
    ofstream out(outfile.c_str());
    if (out) {
      for (unsigned i=R.size()-1; i>=1; i--)
	out << setw(18) << R[i]
	    << setw(18) << -T[i]
	    << setw(18) << -PHI[i]
	    << setw(18) << sqrt(2.0*(E - model->get_pot(R[i])))
	    << setw(18) << model->get_pot(R[i])
	    << endl;
      for (unsigned i=0; i<R.size(); i++)
	out << setw(18) << R[i]
	    << setw(18) << T[i]
	    << setw(18) << PHI[i]
	    << setw(18) <<  sqrt(2.0*(E - model->get_pot(R[i])))
	    << setw(18) << model->get_pot(R[i])
	    << endl;
    }
  }
  
}

void UnboundCoefs::coefs(int L, int M, int Nmax,
			 vector<double>& Times, double dT,
			 AxiSymBiorth *t,
			 vector<CMatrix>& coefs,
			 vector<double>& freqs)
{

// ===================================================================
// Do integrals: begin at Times[min]
// ===================================================================

  Linear1d ts(T, R);
  Linear1d tp(T, PHI);
  Vector pt(1, Nmax);
  CMatrix last(0, freqs.size()-1, 1, Nmax);
  CMatrix curr(0, freqs.size()-1, 1, Nmax);
  KComplex I(0.0, 1.0);

  coefs = vector<CMatrix> (Times.size());
  for (unsigned it=0; it<Times.size();it++) {
    coefs[it].setsize(0, freqs.size()-1, 1, Nmax);
    coefs[it].zero();
  }

  for (unsigned it=0; it<Times.size()-1; it++) {

				// Accumulate
    coefs[it+1] += coefs[it];

				// First integrand
    double Time = Times[it];
    double rr = ts.eval(Time);
    double pp = tp.eval(Time);
    CVector cpt;

    t->potl(Nmax, L, rr, pt);
    cpt = pt;
    for (unsigned iv=0; iv<freqs.size(); iv++) {
      curr[iv] = exp(I*(freqs[iv]*Time + pp*M))*cpt;
    }

				// Trapezoidal rule
    while (Time + dT < Times[it+1]) {
      last = curr;
      
      Time += dT;
      rr = ts.eval(Time);
      pp = tp.eval(Time);

      t->potl(Nmax, L, rr, pt);
      cpt = pt;
      for (unsigned iv=0; iv<freqs.size(); iv++) {
	curr[iv] = exp(I*(freqs[iv]*Time + pp*M))*cpt;
      }

      coefs[it+1] += dT*0.5*(curr + last);
    }

    
				// Last integrand
    last = curr;
    double lastT = Time;

    Time = Times[it+1];
    rr = ts.eval(Time);
    pp = tp.eval(Time);

    t->potl(Nmax, L, rr, pt);
    cpt = pt;
    for (unsigned iv=0; iv<freqs.size(); iv++) {
      curr[iv] = exp(I*(freqs[iv]*Time + pp*M))*cpt;
    }

    coefs[it+1] += (Time-lastT)*0.5*(curr + last);
  }

}
