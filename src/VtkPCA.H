#ifndef _VtkPCA_H
#define _VtkPCA_H

#include <boost/shared_ptr.hpp>
#include <Vector.h>

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
  std::vector<vtkFloatArrayP> vals, vecs;
  // Labels
  std::vector<std::string>    elab;
  // Dataset
  vtkRectilinearGridP         dataSet;
  // Save smoothed basis
  bool                        smooth;
public:
  //! Constructor
  VtkPCA(int N, bool smooth=true);

  //! Add data for two-dimensional cylindrical basis
  void Add(const Vector& eval, const Matrix& evec, int m);

  //! Add data for three-dimensional spherical basis
  void Add(const Vector& eval, const Matrix& evec, int l, int m,
	   char tag);
  
  //! Write output file
  void Write(const std::string& name);
};

typedef boost::shared_ptr<VtkPCA> VtkPCAptr;

#endif
