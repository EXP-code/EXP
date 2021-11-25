#ifndef _InteractSelect_H
#define _InteractSelect_H

#include <random>
#include <memory>

#include <Cspline.H>
#include <Ion.H>

using namespace std;

typedef std::uniform_real_distribution<> dist_type;
typedef std::minstd_rand gen_type;
typedef std::shared_ptr<gen_type> gen_ptr;
typedef std::shared_ptr<dist_type> dist_ptr;

class InteractSelect 
{
private:
  void setSeed();
  gen_type generator;
  dist_ptr uni_dist;
  gen_ptr uni;

public:
  //! Constructor
  InteractSelect();

  //! Return cross section and excitation energy, randomly selected
  //! from all possible excitations by total cross section
  std::pair<double, double> selectCEInteract
  (const std::shared_ptr<Ion> a, const Ion::collType& cumCross);

  //! Compute the free-free radiated energy based on the electron velocity
  double selectFFInteract(const std::pair<double, double>& d) { return d.second; }

  //! Return the lost ionization energy and adjust ionization state
  double DIInterLoss(const std::shared_ptr<Ion> a);

  //! 
  double selectRRInteract(const std::shared_ptr<Ion> a,
			  const std::vector<double>& cumCross, double Ee);
};

#endif
