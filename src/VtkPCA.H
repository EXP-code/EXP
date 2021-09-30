#ifndef _VtkPCA_H
#define _VtkPCA_H

#include <sstream>
#include <memory>
#include <Eigen/Eigen>

//
// VTK stuff
//
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkStructuredPoints.h>
#include <vtkRectilinearGrid.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkLookupTable.h>
#include <vtkVersion.h>

using vtkRectilinearGridP       = vtkSmartPointer<vtkRectilinearGrid>;
using vtkRectilinearGridWriterP = vtkSmartPointer<vtkXMLRectilinearGridWriter>;
using vtkFloatArrayP            = vtkSmartPointer<vtkFloatArray>;

class VtkPCA
{
private:
  int nmax;

  // Knots
  vtkFloatArrayP XX, YY, ZZ;

  // Fields
  std::vector<vtkFloatArrayP> coef, hall, eval, snrv, vecs, vec2, covr;

  // Labels
  std::vector<std::string>    elab;

  // Dataset
  vtkRectilinearGridP         dataSet;

  // Save smoothed basis
  bool                        smooth;

  // Rorder by S/N
  bool                        reorder;
  
public:
  //! Constructor
  VtkPCA(int N, bool reorder=false, bool smooth=false);

  //! Add data for two-dimensional cylindrical basis
  void Add(const Eigen::VectorXd& coef,
	   const Eigen::VectorXd& hall,
	   const Eigen::VectorXd& snr,
	   const Eigen::VectorXd& eval,
	   const Eigen::MatrixXd& evec,
	   const Eigen::MatrixXd& covr,
	   int m);

  //! Add data for three-dimensional spherical basis
  void Add(const Eigen::VectorXd& coef,
	   const Eigen::VectorXd& hall,
	   const Eigen::VectorXd& snr,
	   const Eigen::VectorXd& eval,
	   const Eigen::MatrixXd& evec,
	   const Eigen::MatrixXd& covr,
	   int l, int m);
  
  //! Write output file
  void Write(const std::string& name);
};

typedef std::shared_ptr<VtkPCA> VtkPCAptr;

#endif
