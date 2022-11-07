#ifndef _P2QUANTILE_H_
#define _P2QUANTILE_H_

#include <array>

/**
   Implementation of the quantile estimator from:

   Jain, Raj, and Imrich Chlamtac. “The P² algorithm for dynamic
   calculation of quantiles and histograms without storing
   observations.” Communications of the ACM 28, no. 10 (1985):
   1076-1085.
 */
class P2Quantile
{
private:

  //! Probability
  double p;

  //! Marker indices
  std::array<int, 5> n;

  //! Desired marker positions
  std::array<double, 5> ns;	

  //! Marker updates
  std::array<double, 5> dns;

  //! Marker heights
  std::array<double, 5> q;

  //! Update count so far
  size_t count;

  double Parabolic(int i, double d)
    {
      return q[i] + d / (n[i + 1] - n[i - 1]) *
	( (n[i] - n[i - 1] + d) * (q[i + 1] - q[i]) / (n[i + 1] - n[i]) +
	  (n[i + 1] - n[i] - d) * (q[i] - q[i - 1]) / (n[i] - n[i - 1])
	  );
    }

  double Linear(int i, int d)
  {
    return q[i] + d * (q[i + d] - q[i]) / (n[i + d] - n[i]);
  }

public:

  //! Quantile constructor where @param prob is the quantile
  P2Quantile(double prob=0.5) : p(prob), count(0) {}

  //! Change the probability value and reset
  void reset(double prob)
  {
    p = prob;
    count = 0;
  }

  //! Add a new value
  void addValue(double x)
  {
    if (count < 5) {
      q[count++] = x;
      if (count == 5) {
	std::sort(q.begin(), q.end());
	    
	for (int i=0; i<5; i++) n[i] = i;
	    
	ns[0] = 0;
	ns[1] = 2 * p;
	ns[2] = 4 * p;
	ns[3] = 2 + 2 * p;
	ns[4] = 4;
	    
	dns[0] = 0;
	dns[1] = p / 2;
	dns[2] = p;
	dns[3] = (1 + p) / 2;
	dns[4] = 1;
      }
      
      return;
    }

    int k = 0;
    if (x < q[0]) { q[0] = x; k = 0; }
    else if (x < q[1]) k = 0;
    else if (x < q[2]) k = 1;
    else if (x < q[3]) k = 2;
    else if (x < q[4]) k = 3;
    else { q[4] = x; k = 3; }
    
    for (int i=k+1; i<5; i++) n[i]++;
    for (int i=0; i<5; i++) ns[i] += dns[i];
    
    for (int i=1; i<=3; i++) {
      double d = ns[i] - n[i];
      if (d>=1 && n[i+1]-n[i]>1 || d<=-1 && n[i-1] - n[i]<-1) {
	int dInt = 1;
	if (d==0) dInt =  0;
	if (d <0) dInt = -1;
	double qs = Parabolic(i, dInt);
	if (q[i-1] < qs && qs < q[i+1])
	  q[i] = qs;
	else
	  q[i] = Linear(i, dInt);
	n[i] += dInt;
      }
    }
    
    count++;
  }
    
  //! Get the current quantile value
  double getQuantile()
  {
    if (count == 0)
      throw std::runtime_error("Sequence contains no elements");

    if (count <= 5) {
      std::sort(q.begin(), q.begin()+count);
      int index = std::round((count - 1) * p);
      return q[index];
    }
    
    return q[2];
  }

private:

};

#endif
