// This may look like C code, but it is really -*- C++ -*-

#include <vector>
#include <deque>

#include <Vector.h>
#include <CVector.h>

/**
   Compute time series of coefficients for an arbitrary orbit
 */
class TimeSeriesCoefs
{
private:
  deque<double> R, T, PHI;

  double E, Rperi;
  double delT, Tmax, lnL, M;
  string outfile;

  vector<double> compute_df(double mass, double logL, 
			    vector<double>& ps, AxiSymModel* model);

public:

  //! Constructor
  TimeSeriesCoefs (double Energy, double rperi, double rsoft,
		   double dT, double T, AxiSymModel *model, 
		   double mass, double logL, string OUTFILE="");
  
  //! Do the calculation
  void coefs(int L, int M, int Nmax, int NINT, AxiSymBiorth *t, 
	     vector<KComplex>& Freqs,
	     vector<double>& Times, 
	     vector<CMatrix>& cofs);

};

