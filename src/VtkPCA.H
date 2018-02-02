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

public:
  //! Constructor
  VtkPCA(int N);

  //! Add data
  void AddMatrix(const Vector& eval, const Matrix& evec, int indx);

  //! Write output file
  void Write(const std::string& name);
};

typedef boost::shared_ptr<VtkPCA> VtkPCAptr;

#endif
