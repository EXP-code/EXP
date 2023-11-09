#ifndef _EULER_H
#define _EULER_H

#include <Eigen/Eigen>

Eigen::Matrix3d
return_euler_slater(double PHI, double THETA, double PSI, int BODY);

Eigen::Matrix3d
return_euler(double PHI, double THETA, double PSI, int BODY);


#endif
