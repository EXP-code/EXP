#ifndef _Field_Generator_H
#define _Field_Generator_H

#include <vector>
#include <map>

#include <BasisFactory.H>
#include <Coefficients.H>

namespace Field
{
  
  //! Compute 2d slices and 3d volumes for a collection of
  //! coefficients and a given basis type.
  class FieldGenerator
  {
    
  private:
    
    std::vector<double> times, pmin, pmax;
    std::vector<int>    grid;
    
    //! Sanity check time vector with coefficient DB
    void check_times(Coefs::CoefsPtr coefs);
    
    //! Using MPI
    bool use_mpi;

  public:
    
    //! Constructor
    FieldGenerator(const std::vector<double> &time,
		   const std::vector<double> &pmin,
		   const std::vector<double> &pmax,
		   const std::vector<int>    &grid);

    /** Get a field slices as a map in time and type

	For example:
	.
	.
	// Generate the fields for all coefficients in 'coefs'
	auto db = slices(basis, coefs);

	// Get an image at Time=3.14 for density ("d")
	Eigen::MatrixXf surface = db["3.14"]["d"];
    */
    std::map<double, std::map<std::string, Eigen::MatrixXf>>
    slices(Basis::BasisPtr basis, Coefs::CoefsPtr coefs);
    
    //! Compute surface histogram from particles
    std::map<std::string, Eigen::MatrixXf>
    histogram(PR::PRptr reader, double length, int bins,
	      std::vector<double> center={0.0, 0.0, 0.0});

    //! Write field slices to files.  This will be VTK your build is
    //! compiled with VTK and ascii tables otherwise.
    void file_slices(Basis::BasisPtr basis, Coefs::CoefsPtr coefs,
		    const std::string prefix, const std::string outdir=".");

    /** Get a field volumes as a map in time and type.
    
	Usage as for slices but returns an Eigen::Tensor rather than an 
	Eigen::MatrixXf.
    */
    std::map<double, std::map<std::string, Eigen::Tensor<float, 3>>>
    volumes(Basis::BasisPtr basis, Coefs::CoefsPtr coefs);
    
    //! Write field slices to files.  This will be VTK your build is
    //! compiled with VTK and ascii tables otherwise.
    void file_volumes(Basis::BasisPtr basis, Coefs::CoefsPtr coefs,
		      const std::string prefix, const std::string outdir=".");
    //@}
    
  };

}
// END namespace Field
  
#endif

  
