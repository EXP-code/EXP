#ifndef _PseudoAccel_H_
#define _PseudoAccel_H_

#include <deque>
#include <array>

#include <Eigen/Eigen>

#include <QuadLS.H>

class PseudoAccel
{
private:

  //! Maximum number of time points
  unsigned int Nsize;

  //! Type flags
  bool CENTER=false, AXIS=false;

  //! Queue of elements
  using Elem = std::array<double, 7>;
  std::deque<Elem> queue;

  //! Storage 
  Eigen::Vector3d accel = Eigen::Vector3d::Zero();
  Eigen::Vector3d omega = Eigen::Vector3d::Zero();
  Eigen::Vector3d domdt = Eigen::Vector3d::Zero();

public:

  //! Construct a pseudo acceleration estimator
  PseudoAccel(unsigned int Nsize, bool center, bool axis) :
    Nsize(Nsize), CENTER(center), AXIS(axis) {}

  //! Add a center element to the pseudo accelration estimator
  void add(double t, const Eigen::Vector3d& c, const Eigen::Vector3d& a) {
    queue.push_back({t, c(0), c(1), c(2), a(0), a(1), a(2)});
    if (queue.size() > Nsize) queue.pop_front();
  }
	  
  std::tuple<Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector3d&> operator()()
  {
    if (CENTER or AXIS) {

      std::vector<double> t, x, y, z, u, v, w;
      for (auto &e : queue) {
	t.push_back(e[0]);
	if (CENTER) {
	  x.push_back(e[1]);
	  y.push_back(e[2]);
	  z.push_back(e[3]);
	}
	if (AXIS) {
	  u.push_back(e[4]);
	  v.push_back(e[5]);
	  w.push_back(e[6]);
	}
      }

      if (CENTER) {
	accel << 
	  2.0*std::get<0>(QuadLS(t, x).coefs()),
	  2.0*std::get<0>(QuadLS(t, y).coefs()),
	  2.0*std::get<0>(QuadLS(t, z).coefs());
      }

      if (AXIS) {
	auto [a, b, c] = QuadLS(t, u).coefs();
	auto [d, e, f] = QuadLS(t, v).coefs();
	auto [g, h, i] = QuadLS(t, w).coefs();
	
	double T = t.back();	// Last eval time

	Eigen::Vector3d n, dndt, d2ndt2;
	n      << a*T*T+b*T+c, d*T*T+e*T+f, g*T*T+h*T+i;
	dndt   << 2.0*a*T+b, 2.0*d*T+e, 2.0*g*T+h;
	d2ndt2 << 2.0*a, 2.0*d, 2.0*g;

	// Get angular acceleration and its derivative
	omega = n.cross(dndt);
	domdt = n.cross(d2ndt2);
      }

    }

    return {accel, omega, domdt};
  }
};

#endif
