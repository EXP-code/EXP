#include <iostream>
#include <iomanip>
#include <Eigen/Eigen>


int main()
{
  Eigen::MatrixXd M(3, 3);
  Eigen::VectorXd e(3);
  
  M(0, 0) = 1.0;
  M(0, 1) = 0.5;
  M(0, 2) = 0.1;
  
  M(1, 1) = 1.0;
  M(1, 2) = 0.2;
  
  M(2, 2) = 1.0;

  M(1, 0) = M(0, 1);
  M(2, 0) = M(0, 2);
  M(2, 1) = M(1, 2);
  
  EigenSolver<MatrixXd> es(M, false);
  std::cout << "The eigenvalues of the 3x3 matrix of ones are:" 
	    << std::endl << es.eigenvalues() << std::endl;  
}
	



