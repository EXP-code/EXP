#ifndef _QUADLS_H_
#define _QUADLS_H_

#include <iostream>
#include <iomanip>
#include <tuple>

// Take a pair of data arrays (x, y) and perform a quadratic regression
template <typename T>
class QuadLS
{
private:
  double _a, _b, _c;

  void fit( const T &x, const T &y )
  {
    auto n = x.size();

    if (n==y.size() and n>0) {

      double sumx=0, sumy=0, sumxy=0, sumx2y=0, sumx2=0, sumx3=0, sumx4=0;

      for (int i=0; i<n; i++) {
	sumx   += x[i];
	sumy   += y[i];
	sumx2  += x[i]*x[i];
	sumxy  += x[i]*y[i];
	sumx2y += x[i]*x[i]*y[i];
	sumx3  += x[i]*x[i]*x[i];
	sumx4  += x[i]*x[i]*x[i]*x[i];
      }
      
      double Sxx   = sumx2  - sumx*sumx/n;
      double Sxy   = sumxy  - sumx*sumy/n;
      double Sxx2  = sumx3  - sumx*sumx2/n;
      double Sx2y  = sumx2y - sumx2*sumy/n;
      double Sx2x2 = sumx4  - sumx2*sumx2/n;
      
      double denom = Sxx*Sx2x2 - Sxx2*Sxx2;

      if (fabs(denom)>0.0) {
	_a = (Sx2y*Sxx  - Sxy*Sxx2 ) / denom;
	_b = (Sxy*Sx2x2 - Sx2y*Sxx2) / denom;
	_c = (sumy - sumx2*_a - sumx*_b) / n;
      } else _a = _b = _c = 0.0;

    }
    else {
      _a = _b = _c = 0.0;
    }
  }

public:

  //! Null constructor
  QuadLS() { }

  //! Constructor
  QuadLS( const T &x, const T &y )
  { fit(x, y); }
	
  //! Destructor
  virtual ~QuadLS() { }

  //! Access coefficients
  std::tuple<double, double, double> coefs() const
  { return std::make_tuple( _a, _b, _c ); }

  //! Evaluate the quadratic function
  double operator()( const double x ) const
  { return _a * x * x + _b * x + _c; }
};

#endif
