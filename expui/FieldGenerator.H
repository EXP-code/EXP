#ifndef _Field_Generator_H
#define _Field_Generator_H

#include <string_view>
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
    Eigen::MatrixXd     mesh;
    
    //! Sanity check time vector with coefficient DB
    void check_times(CoefClasses::CoefsPtr coefs);
    
    //! Using MPI
    bool use_mpi = false;

    //! Perform midplane evaluation
    bool midplane = false;

    //! Midplane search height
    double colheight = 4.0;

  public:
    
    //! Null constructor
    FieldGenerator() = default;

    //! Constructor for a rectangular grid
    FieldGenerator(const std::vector<double> &time,
		   const std::vector<double> &pmin,
		   const std::vector<double> &pmax,
		   const std::vector<int>    &grid);

    //! Constructor for an arbitrary point mesh.  The mesh should be
    //! an Nx3 array
    FieldGenerator(const std::vector<double> &time,
		   const Eigen::MatrixXd     &mesh);

    /** Get field quantities at user defined points

	For example:
	.
	.
	// Generate the fields for all coefficients in 'coefs'
	auto db = points(basis, coefs);

	// Get fields evaluated at Time=3.14 and for all points in
	// your supplied mesh for density ("dens")
	Eigen::MatrixXf points = db["3.14"]["dens"];
    */
    std::map<double, std::map<std::string, Eigen::VectorXf>>
    points(BasisClasses::BasisPtr basis, CoefClasses::CoefsPtr coefs);

    /** Get a field slices as a map in time and type

	For example:
	.
	.
	// Generate the fields for all coefficients in 'coefs'
	auto db = slices(basis, coefs);

	// Get an image at Time=3.14 for density ("dens")
	Eigen::MatrixXf surface = db["3.14"]["dens"];
    */
    std::map<double, std::map<std::string, Eigen::MatrixXf>>
    slices(BasisClasses::BasisPtr basis, CoefClasses::CoefsPtr coefs);
    
    //@{
    /** Get a field probes as a map in time and type

	For example:
	.
	.
	// Generate the fields for all coefficients in 'coefs' along the 
	// line specfied by (beg, end) where beg, end are three vectors 
	// describing points in space and num is the number of evaluations
	auto db = lines(basis, coefs, beg, end, num);

	// The code will provide three additional vectors: x, y, z,
	// arc where x, y, z, are the evaulation positions along the
	// line and arc describing the distance along the line from
	// beg to end

	// Get an image at Time=3.14 for density ("dens")
	arc = db["3.14"]["arc"];
	pot = db["3.14"]["potl"] - db["3.14"]["potl m=0"];
    */
    std::map<double, std::map<std::string, Eigen::VectorXf>>
    lines(BasisClasses::BasisPtr basis, CoefClasses::CoefsPtr coefs,
	  std::vector<double> beg, std::vector<double> end, int num);
    
    void
    file_lines(BasisClasses::BasisPtr basis, CoefClasses::CoefsPtr coefs,
	       std::vector<double> beg, std::vector<double> end, int num,
	       const std::string prefix, const std::string outdir=".");    
    //@}

    //! Compute surface histogram from particles
    std::map<std::string, Eigen::MatrixXf>
    histogram2d(PR::PRptr reader, std::vector<double> center={0.0, 0.0, 0.0});

    //! Compute histogram from particles
    Eigen::VectorXf
    histogram1d(PR::PRptr reader, double rmax, int nbins, std::string proj,
		std::vector<double> center={0.0, 0.0, 0.0});

    //! Compute spherical log histogram from particles
    std::tuple<Eigen::VectorXf, Eigen::VectorXf>
    histo1dlog(PR::PRptr reader, double rmin, double rmax, int nbins,
		std::vector<double> center={0.0, 0.0, 0.0});

    //! Write field slices to files.  This will be VTK your build is
    //! compiled with VTK and ascii tables otherwise.
    void file_slices(BasisClasses::BasisPtr basis, CoefClasses::CoefsPtr coefs,
		     const std::string prefix, const std::string outdir=".");

    /** Get a field volumes as a map in time and type.
    
	Usage as for slices but returns an Eigen::Tensor rather than an 
	Eigen::MatrixXf.
    */
    std::map<double, std::map<std::string, Eigen::Tensor<float, 3>>>
    volumes(BasisClasses::BasisPtr basis, CoefClasses::CoefsPtr coefs);
    
    //! Write field slices to files.  This will be VTK your build is
    //! compiled with VTK and ascii tables otherwise.
    void file_volumes(BasisClasses::BasisPtr basis, CoefClasses::CoefsPtr coefs,
		      const std::string prefix, const std::string outdir=".");
    //@}
    
    //! Turn on/off midplane evaluation (only effective for disk basis
    //! and slices)
    void setMidplane(bool value) { midplane = value; }

    //! Height above/below the plane for midplane search in disk scale
    //! lengths
    void setColumnHeight(double value) { colheight = value; }

  };

}
// END namespace Field
  
#endif

  
