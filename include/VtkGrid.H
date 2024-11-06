#ifndef _VtkGrid_H
#define _VtkGrid_H

#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <memory>
#include <limits>
#include <map>

#ifdef HAVE_VTK

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

#include <ThreeDGrid.H>

/**
   A ThreeDGrid implementation of the rectangular grid database
   implement in VTK and designed for consumption by PyVTK or ParaView.
 */
class VtkGrid : public ThreeDGrid
{
private:
  // Knots
  vtkSmartPointer<vtkFloatArrayP> XX, YY, ZZ;

  // Dataset
  vtkRectilinearGridP dataSet;

public:
  //! Constructor
  VtkGrid(int nx, int ny, int nz,
	 double xmin, double xmax,
	 double ymin, double ymax,
	 double zmin, double zmax);

  //! Add data for two-dimensional cylindrical basis
  void Add(const std::vector<double>& data, const std::string& name);

  //! Write output file
  void Write(const std::string& name);
};

#else

#include <ThreeDGrid.H>

/**
   A ThreeDGrid implementation of the rectangular grid database
   implement in VTK and designed for consumption by PyVTK or ParaView.
 */
class VtkGrid : public ThreeDGrid
{
private:
  // Grid coordinates
  std::map<std::string, std::vector<float>> coord;

  // Dataset
  std::map<std::string, std::vector<float>> dataSet;

  // Write header
  void writeBeg(std::ofstream& file);

  // Write footer
  void writeEnd(std::ofstream& file);

  // Write scalar fields
  void writeFields(std::ofstream& file);

  // Write coordinates
  void writeCoords(std::ofstream& file);

  // String replacement utility
  void replaceAll(std::string& str,
		  const std::string& from,
		  const std::string& to);

public:
  //! Constructor
  VtkGrid(int nx, int ny, int nz,
	 double xmin, double xmax,
	 double ymin, double ymax,
	 double zmin, double zmax);

  //! Add data for two-dimensional cylindrical basis
  void Add(const std::vector<double>& data, const std::string& name);

  //! Write output file
  void Write(const std::string& name);
};

#endif

typedef std::shared_ptr<VtkGrid> VtkGridPtr;

#endif

