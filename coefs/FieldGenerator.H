#ifndef _Field_Generator_H
#define _Field_Generator_H

#include <vector>
#include <map>

#include <BasisFactory.H>
#include <CoefFactory.H>

namespace Field
{
  
  class FieldGenerator
  {
    
  private:
    
    std::vector<double> time, pmin, pmax;
    std::vector<int>    grid;
    
    //! Sanity check time vector with coefficient DB
    void check_times(Coefs::CoefsPtr coefs);
    
  public:
    
    //! Constructor
    FieldGenerator(const std::vector<double> &time,
		   const std::vector<double> &pmin,
		   const std::vector<double> &pmax,
		   const std::vector<int>    &grid);
    //@{
    //! Get a field slice
    std::map<double, std::map<std::string, Eigen::MatrixXd>>
    slices(Basis::BasisPtr basis, Coefs::CoefsPtr coefs);
    
    void vtk_slices(Basis::BasisPtr basis, Coefs::CoefsPtr coefs);
    //@}
    
    //@{
    //! Get a field volume
    std::map<double, std::map<std::string, Eigen::Tensor<float, 3>>>
    volumes(Basis::BasisPtr basis, Coefs::CoefsPtr coefs);
    
    void vtk_volumes(Basis::BasisPtr basis, Coefs::CoefsPtr coefs);
    //@}
    
  };

}
// END namespace Field
  
#endif

  
