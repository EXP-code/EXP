#ifndef _UserWake_H
#define _UserWake_H

#include <vector>
#include <string>

/** Print the wake for the desired components

    @param NUMX is the number of points in the X-direction
    @param NUMY is the number of points in the Y-direction
    @param XMIN is the minimum point in the X-direction
    @param XMAX is the maximum point in the X-direction
    @param YMIN is the minimum point in the Y-direction
    @param YMAX is the maximum point in the Y-direction
    @param PHI is the Euler angle
    @param PSI is the Euler angle
    @param THETA is the Euler angle
    @param NSTEP is the number of steps between frames
    @param filename is the suffix of the output file
    @param C(i) is a component to include in the wake.  The index i must begin at 1 and increase by 1 for each additional component.
*/
class UserWake : public ExternalForce
{
private:
  
  vector<string> C;
  vector<Component*> c0;
  vector<double> r, theta, phi;
  vector< vector<float> > data0, data1;

  bool ORIENT, first;

  string filename;

  int NUMX, NUMY;
  double XMIN, XMAX, YMIN, YMAX;
  double PHI, PSI, THETA;

  vector<string> names;
  int NSTEP, numComp, count, nlast, nnext, nbeg, nend, npix;

  void determine_acceleration_and_potential(void);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  void userinfo();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  UserWake(const YAML::Node& conf);

  //! Destructor
  ~UserWake();

};

Eigen::Matrix3d return_euler_slater(double, double, double, int);


#endif
