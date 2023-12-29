// exputils signatures

#ifndef _EXPUTILS_H
#define _EXPUTILS_H

#include <string>
#include <Eigen/Eigen>

//! Compare biorthogoal matrices with idenity and return worst case
//! error and success value based on __EXP__::orthoTol
std::tuple<bool, double, std::vector<double>>
orthoCompute(const std::vector<Eigen::MatrixXd>& tests);

//! Use orthoCompute() to compare biorthogonal matrices with idenity
//! and throw std::runtime_error if worst case exceeds tolerance
void
orthoCompare(const std::vector<Eigen::MatrixXd>& tests,
	     const std::string& classname, const std::string& indexname);

#endif
