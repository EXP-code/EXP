#include <iostream>
#include <fstream>
#include <string>

#include <interp.H>
#include <massmodel.H>
#include <biorth.H>
#include <TimeSeriesCoefs.H>

vector<double> TimeSeriesCoefs::compute_df(double smass, double lnL,
					   vector<double>& ps, AxiSymModPtr m)
{
  double rmin = m->get_min_radius();
  double r = sqrt(ps[0]*ps[0] + ps[1]*ps[1] + rmin*rmin);
  double v = sqrt(ps[2]*ps[2] + ps[3]*ps[3]) + 1.0e-18;

  // From usual Chandra formula:
  //
  // 4*pi*[erf(x) - 2x/sqrt(pi)*exp(-x*x)]
  // assuming that x = v/vc
  //
  double x   = v/sqrt(m->get_mass(r)/r);
  double cof = 4.0*M_PI*(erf(x) - 2.0*x/sqrt(M_PI)*exp(-x*x));

  double accel = -cof*lnL*smass*m->get_density(r)/(v*v);

  vector<double> ret(2, 0.0);
  ret[0] = accel*ps[2]/v;
  ret[1] = accel*ps[3]/v;

  return ret;
}

TimeSeriesCoefs::TimeSeriesCoefs(double Energy, double rperi, double rsoft,
				 double dT, double Time, AxiSymModPtr model, 
				 double M, double lnL, string OUTFILE) :
  E(Energy), Rperi(rperi), delT(dT), Tmax(Time)
{

  double deltaE = E - model->get_pot(Rperi);
  double rmin   = model->get_min_radius();
  if (deltaE<0.0) {
    cerr << "Rperi=" << Rperi << ">Rapo for E=" << E << endl;
    exit(-1);
  }

  // =================================
  // Compute orbit
  // =================================

  double VTperi = sqrt(2.0*deltaE), r, rs, dPot;

  deque< vector<double> > PS;
  vector<double> cur(4), lst(4), frc(4), drg(2, 0.0);
  
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
    
    r    = sqrt(cur[0]*cur[0] + cur[1]*cur[1] + rmin*rmin);
    rs   = r*r + rsoft*rsoft;
    dPot = model->get_mass(r)/rs;

    if (lnL > 0.0) drg = compute_df(M, lnL, cur, model);

    frc[0] = cur[2];
    frc[1] = cur[3];
    frc[2] = -dPot*cur[0]/r + drg[0];
    frc[3] = -dPot*cur[1]/r + drg[1];

    for (int i=0; i<4; i++) lst[i] = cur[i] + 0.5*delT*frc[i];

    r = sqrt(lst[0]*lst[0] + lst[1]*lst[1] + rmin*rmin);
    rs   = r*r + rsoft*rsoft;
    dPot = model->get_mass(r)/rs;

    if (lnL > 0.0) drg = compute_df(M, lnL, lst, model);

    frc[0] = lst[2];
    frc[1] = lst[3];
    frc[2] = -dPot*lst[0]/r + drg[0];
    frc[3] = -dPot*lst[1]/r + drg[1];

    for (int i=0; i<4; i++) cur[i] += delT*frc[i];

    time += delT;

    T.push_back(time);
    R.push_back(sqrt(cur[0]*cur[0] + cur[1]*cur[1]));
    PHI.push_back(atan2(cur[1], cur[0]));
    PS.push_back(cur);
  }
  
  time = 0.0;

  cur[0] = Rperi;
  cur[1] = 0.0;
  cur[2] = 0.0;
  cur[3] = VTperi;

  while (time > -Tmax) {
    
    r    = sqrt(cur[0]*cur[0] + cur[1]*cur[1] + rmin*rmin);
    rs   = r*r + rsoft*rsoft;
    dPot = model->get_mass(r)/rs;

    if (M>0.0) drg = compute_df(M, lnL, cur, model);

    frc[0] = cur[2];
    frc[1] = cur[3];
    frc[2] = -dPot*cur[0]/r + drg[0];
    frc[3] = -dPot*cur[1]/r + drg[1];

    for (int i=0; i<4; i++) lst[i] = cur[i] - 0.5*delT*frc[i];

    r = sqrt(lst[0]*lst[0] + lst[1]*lst[1]) + 1.0e-18;
    dPot = model->get_dpot(r);

    if (M>0.0) drg = compute_df(M, lnL, lst, model);

    frc[0] = lst[2];
    frc[1] = lst[3];
    frc[2] = -dPot*lst[0]/r + drg[0];
    frc[3] = -dPot*lst[1]/r + drg[1];

    for (int i=0; i<4; i++) cur[i] -= delT*frc[i];

    time -= delT;

    T.push_front(time);
    R.push_front(sqrt(cur[0]*cur[0] + cur[1]*cur[1]));
    PHI.push_front(atan2(cur[1], cur[0]));
    PS.push_front(cur);
  }


  if (OUTFILE.size()) {

    std::string file = OUTFILE + ".orbit";
    std::ofstream out(file);
    if (out) {
      for (unsigned i=0; i<PS.size(); i++)
	out << setw(18) << T[i]
	    << setw(18) << R[i]
	    << setw(18) << PHI[i]
	    << setw(18) << PS[i][0]
	    << setw(18) << PS[i][1]
	    << setw(18) << PS[i][2]
	    << setw(18) << PS[i][3]
	    << endl;
    }
  }
}

void TimeSeriesCoefs::coefs(int L, int M, int Nmax, int NINT,
			    AxiSymBioPtr t,
			    std::vector<complex<double>>& Freqs,
			    std::vector<double>& Times,
			    std::vector<Eigen::MatrixXcd>& coefs)
{
  constexpr std::complex<double> I(0.0, 1.0);

// ===================================================================
// Do integrals: begin at Times[min]
// ===================================================================

  coefs.resize(Times.size());
  for (auto & v : coefs) {
    v.resize(Freqs.size(), Nmax);
    v.setZero();
  }

  LegeQuad lq(NINT);
  Eigen::VectorXd pt(Nmax);
  Eigen::VectorXcd cpt;

  double rr, pp, Tmin=Times.front(), Tmax=Times.back();

  for (unsigned it=0; it<Times.size(); it++) {
    double tt = Times[it];
    for (unsigned iv=0; iv<Freqs.size(); iv++) {
      double Time = Tmin;
      while (Time<Tmax) {
	double dT = min<double>(2.0*M_PI/std::abs(Freqs[iv]), Tmax-Time);
	for (int jt=0; jt<NINT; jt++) {
	  double time = Time + dT*lq.knot(jt);
	  rr = odd2(time, T, R);
	  pp = odd2(time, T, PHI);
	  t->potl(Nmax, L, rr, pt);
	  cpt = pt;
	  coefs[it].row(iv) += cpt * exp(I*Freqs[iv]*(tt-time)-I*pp*static_cast<double>(M)) * 
	    dT*lq.weight(jt);
	}
	Time += dT;
      }
    }
  }
}
