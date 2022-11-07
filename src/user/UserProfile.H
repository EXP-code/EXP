#ifndef _UserWake_H
#define _UserWake_H

#include <vector>
#include <string>

/** Print the spherical component of the profile for a particular component

    @param NUMR is the number of radial points
    @param RMIN is the minimum radius
    @param RMAX is the maximum radius
    @param LOGR true means logarithmic spacing between RMIN and RMAX
    @param NTHETA is the number of colatitude knots
    @param NPHI is the number of azimuthal knots
    @param DT is the time interval between evaluations
    @param filename is the suffix of the output file
    @param C(i) is a component to include in the wake.  The index i must begin at 1 and increase by 1 for each additional component.
*/
class UserProfile : public ExternalForce
{
private:
  
  vector<string> C;
  vector<Component*> c0;
  vector<double> rho, mass, pot, rho1, mass1, pot1;

  bool LOGR, first;

  string filename;

  int NUMR, NTHETA, NPHI;
  double RMIN, RMAX, dR;

  double DT, tnext;
  int numComp, count;
  vector<int> nbeg, nend, numb;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg) { return 0x0; }
  void initialize();

  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserProfile(const YAML::Node& conf);

  //! Destructor
  ~UserProfile();

};

#endif
