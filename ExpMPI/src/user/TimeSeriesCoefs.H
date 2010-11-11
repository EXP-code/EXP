// This may look like C code, but it is really -*- C++ -*-

#include <Vector.h>
#include <CVector.h>

/**
   Compute time series of coefficients for an arbitrary orbit
 */
class TimeSeriesCoefs
{
private:
  vector<double> R, T, PHI;

  double E, Rperi;
  double delT, Tmax;
  string outfile;

public:

  //! Constructor
  TimeSeriesCoefs (double Energy, double rperi, double dT, double T,
		   AxiSymModel *model, string OUTFILE="");
  
  //! Do the calculation
  void coefs(int L, int M, int Nmax, int NINT, AxiSymBiorth *t, 
	     vector<KComplex>& Freqs,
	     vector<double>& Times, 
	     vector<CMatrix>& cofs);

};

