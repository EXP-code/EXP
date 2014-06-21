#ifndef _InteractSelect_H
#define _InteractSelect_H

#include <boost/random.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "Cspline.H"
#include "Ion.H"

using namespace std;

#define planck 6.6260693e-27
#define light 29979245800.
#define boltzmann 1.3806504e-16
#define boltzmannEv 8.617 3324e−5

typedef boost::uniform_real<> dist_type;
typedef boost::variate_generator<boost::minstd_rand, dist_type> gen_type;
typedef boost::shared_ptr<gen_type> gen_ptr;
typedef boost::shared_ptr<dist_type> dist_ptr;

class InteractSelect 
{
private:
  void setSeed();
  boost::minstd_rand generator;
  dist_ptr uni_dist;
  gen_ptr uni;

public:
  //! Constructor
  InteractSelect();

  //! Return excitation energy, randomly selected from all possible
  //! excitations by cross section
  double selectCEInteract
  (const Ion& a, const Ion::collType& cumCross);

  //! Compute the free-free radiated energy based on the electron velocity
  double selectFFInteract(const Ion& a, double E);

  //! Return the lost ionization energy and adjust ionization state
  double DIInterLoss(const chdata& ch, const Ion& a);

  //! 
  double selectRRInteract(const Ion& a, const std::vector<double>& cumCross, double Ee);
};

#endif
