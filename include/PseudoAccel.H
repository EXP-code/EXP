#ifndef _PseudoAccel_H_
#define _PseudoAccel_H_

#include <deque>
#include <array>

#include <Eigen/Eigen>

#include <QuadLS.H>

class PseudoAccel
{
private:

  //! Maximum number of centers
  unsigned int Nsize;

  //! Queue of center elements
  using Elem = std::array<double, 4>;
  std::deque<Elem> queue;

public:

  //! Construct a pseudo acceleration estimator
  PseudoAccel(unsigned int Nsize) : Nsize(Nsize) {}

  //! Add a center element to the pseudo accelration estimator
  void add(double t, const Eigen::Vector3d& c) {
    queue.push_back({t, c(0), c(1), c(2)});
    if (queue.size() > Nsize) queue.pop_front();
  }
	  
  Eigen::Vector3d accel()
  {
    std::vector<double> t, x, y, z;
    for (auto &e : queue) {
      t.push_back(e[0]);
      x.push_back(e[1]);
      y.push_back(e[2]);
      z.push_back(e[3]);
    }

    return {2.0*std::get<0>(QuadLS(t, x).coefs()),
	    2.0*std::get<0>(QuadLS(t, y).coefs()),
	    2.0*std::get<0>(QuadLS(t, z).coefs())};
  }
};

#endif
