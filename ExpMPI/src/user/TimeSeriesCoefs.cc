#include <interp.h>
#include <CVector.h>
// #include <CMatrix.h>
#include <massmodel.h>
#include <biorth.h>
#include <TimeSeriesCoefs.H>

TimeSeriesCoefs::TimeSeriesCoefs(double Energy, double rperi, 
				 double dT, double Time,
				 AxiSymModel *model, string OUTFILE) :
  E(Energy), Rperi(rperi), delT(dT), Tmax(Time)
{

  double deltaE = E - model->get_pot(Rperi);
  if (deltaE<0.0) {
    cerr << "Rperi=" << Rperi << ">Rapo for E=" << E << endl;
    exit(-1);
  }

  // =================================
  // Compute orbit
  // =================================

  double VTperi = sqrt(2.0*deltaE), r, dPot;
  vector< vector<double> > PS;
  vector<double> cur(4), lst(4), frc(4);
  
  cur[0] = Rperi;
  cur[1] = 0.0;
  cur[2] = 0.0;
  cur[3] = VTperi;

  T.push_back(0.0);
  PS.push_back(cur);
  R.push_back(Rperi);
  PHI.push_back(0.0);

  //
  // Midpoint method (aka RK2)
  //

  double time = 0.0;

  while (time < Tmax) {
    
    r = sqrt(cur[0]*cur[0] + cur[1]*cur[1]);
    dPot = model->get_dpot(r);

    frc[0] = cur[2];
    frc[1] = cur[3];
    frc[2] = -dPot*cur[0]/r;
    frc[3] = -dPot*cur[1]/r;

    for (int i=0; i<4; i++) lst[i] = cur[i] + 0.5*delT*frc[i];

    r = sqrt(lst[0]*lst[0] + lst[1]*lst[1]);
    dPot = model->get_dpot(r);

    frc[0] = lst[2];
    frc[1] = lst[3];
    frc[2] = -dPot*lst[0]/r;
    frc[3] = -dPot*lst[1]/r;

    for (int i=0; i<4; i++) cur[i] += delT*frc[i];

    time += delT;

    T.push_back(time);
    R.push_back(sqrt(cur[0]*cur[0] + cur[1]*cur[1]));
    PHI.push_back(atan2(cur[1], cur[0]));
    PS.push_back(cur);
  }
  
  string file = OUTFILE + ".orbit";
  ofstream out(file.c_str());
  if (out) {
    for (unsigned i=PS.size()-1; i>=1; i--)
      out << setw(18) << -T[i]
	  << setw(18) <<  PS[i][0]
	  << setw(18) << -PS[i][1]
	  << setw(18) << -PS[i][2]
	  << setw(18) <<  PS[i][3]
	  << endl;
    for (unsigned i=0; i<PS.size(); i++)
      out << setw(18) << T[i]
	  << setw(18) << PS[i][0]
	  << setw(18) << PS[i][1]
	  << setw(18) << PS[i][2]
	  << setw(18) << PS[i][3]
	  << endl;
  }
}

void TimeSeriesCoefs::coefs(int L, int M, int Nmax, int NINT,
			    AxiSymBiorth *t,
			    vector<KComplex>& Freqs,
			    vector<double>& Times,
			    vector<CMatrix>& coefs)
{
  KComplex I(0.0, 1.0);

// ===================================================================
// Do integrals: begin at Times[min]
// ===================================================================

  vector<double> TT, RR, PP;
  for (unsigned i=R.size()-1; i>=1; i--) {
    TT.push_back(-T[i]);
    RR.push_back(R[i]);
    PP.push_back(-PHI[i]);
  }
  for (unsigned i=0; i<R.size(); i++) {
    TT.push_back(T[i]);
    RR.push_back(R[i]);
    PP.push_back(PHI[i]);
  }
  
  coefs = vector<CMatrix>(Times.size());
  for (unsigned it=0; it<Times.size(); it++) {
    coefs[it] = CMatrix(0, Freqs.size()-1, 1, Nmax);
    coefs[it].zero();
  }

  Vector pt(1, Nmax);
  CVector cpt;
  double rr, pp;
  LegeQuad lq(NINT);

  for (unsigned it=1; it<Times.size(); it++) {
    double tt = Times[it];
    for (unsigned iv=0; iv<Freqs.size(); iv++) {
      double T = Times[0], TT;
      while (T<Times[it]) {
	double dT = min<double>(2.0*M_PI/fabs(Freqs[iv]), tt-T);
	for (int jt=1; jt<=NINT; jt++) {
	TT = T + dT*lq.knot(jt);
	rr = odd2(TT, Times, RR);
	pp = odd2(TT, Times, PP);
	t->potl(Nmax, L, rr, pt);
	cpt = pt;
	coefs[it][iv] += cpt * exp(I*Freqs[iv]*(TT-tt)-I*pp*M) * 
	  dT*lq.weight(jt);
	}
	T += dT;
      }
    }
  }
}
