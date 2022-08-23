#ifndef _Field_Generator_H
#define _Field_Generator_H

#include <vector>
#include <map>

#include <BasisFactory.H>
#include <CoefFactory.H>

class FieldGenerator
{

private:

  std::vector<double> time, pmin, pmax;
  std::vector<int>    grid;

  //! Sanity check time vector with coefficient DB
  void check_times(CoefFactory& coefs);

public:

  //! Constructor
  FieldGenerator(const std::vector<double> &time,
		 const std::vector<double> &pmin,
		 const std::vector<double> &pmax,
		 const std::vector<int>    &grid);
  //@{
  //! Get a field slice
  std::map<double, std::map<std::string, Eigen::MatrixXd>>
  slices(BasisFactory& basis, CoefFactory& coefs);

  void vtk_slices(BasisFactory& basis, CoefFactory& coefs);
  //@}

  //@{
  //! Get a field volume
  std::map<double, std::map<std::string, Eigen::Tensor<float, 3>>>
  volumes(BasisFactory& basis, CoefFactory& coefs);

  void vtk_volumes(BasisFactory& basis, CoefFactory& coefs);
  //@}

};

#endif
