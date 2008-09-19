// This may look like C code, but it is really -*- C++ -*-

#include <Timer.h>
#include <Vector.h>
#include <CVector.h>


/**
   Compute Fourier coefficients for an unbound orbit
 */
class UnboundCoefs
{
private:
  vector<double> R, T, PHI;

  int m;
  int lmax;
  Timer timer;
  double tmax;
  double E, Rperi, Redge, deltaR;

public:

  //! Constructor
  UnboundCoefs (double Energy, double rperi, double redge, double delR,
		AxiSymModel *model, string OUTFILE="");
  
  //! Do the calculation
  void coefs(int L, int M, int Nmax,
	     vector<double>& Times, double dT,
	     AxiSymBiorth *t, vector<CMatrix>& cofs, vector<double>& freqs);

  double Tmax() { return tmax; }

};

