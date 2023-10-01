// exputils signatures

#ifndef _EXPUTILS_H
#define _EXPUTILS_H

#include <string>
#include <Eigen/Eigen>

void orthoTest(const std::vector<Eigen::MatrixXd>& tests,
	       const std::string& classname, const std::string& indexname);

#endif
