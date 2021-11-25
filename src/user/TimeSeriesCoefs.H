#include <vector>
#include <deque>
#include <complex>

#include <Eigen/Eigen>

/**
   Compute time series of coefficients for an arbitrary orbit
 */
class TimeSeriesCoefs
{
private:
  deque<double> R, T, PHI;

  double E, Rperi;
  double delT, Tmax, lnL, M;
  std::string outfile;

  vector<double> compute_df(double mass, double logL, 
			    vector<double>& ps, AxiSymModPtr model);

public:

  //! Constructor
  TimeSeriesCoefs (double Energy, double rperi, double rsoft,
		   double dT, double T, AxiSymModPtr model, 
		   double mass, double logL, string OUTFILE="");
  
  //! Do the calculation
  void coefs(int L, int M, int Nmax, int NINT, AxiSymBioPtr t, 
	     vector<complex<double>>& Freqs,
	     vector<double>& Times, 
	     vector<Eigen::MatrixXcd>& coefs);

};

