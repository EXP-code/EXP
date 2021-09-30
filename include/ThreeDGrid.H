#ifndef _ThreeDGrid_H
#define _ThreeDGrid_H

#include <vector>
#include <sstream>
#include <memory>

/**
   This is an API for populating a 3d grid rectangular grid database
   for diagnostic output.

   The user provides the size of the grid and extents.  Calls to Add()
   add additiional data fields with the predefined rectangular grid
   shape.  The input array is flattened and packed in C (not Fortran)
   order.  A call to Write() dumps the array in the format implemented
   by the derived class.
 */
class ThreeDGrid
{
protected:
  int nx, ny, nz;
  float xmin, xmax, ymin, ymax, zmin, zmax;

public:

  //! Constructor
  ThreeDGrid(int nx, int ny, int nz,
	     double xmin, double xmax,
	     double ymin, double ymax,
	     double zmin, double zmax) :
    nx(nx), ny(ny), nz(nz),
    xmin(xmin), xmax(xmax),
    ymin(ymin), ymax(ymax),
    zmin(zmin), zmax(zmax) {}
  
  //! Null constructor
  ThreeDGrid() : nx(0), ny(0), nz(0),
		 xmin(0), xmax(0),
		 ymin(0), ymax(0),
		 zmin(0), zmax(0) {}

  //! Add a data field to the rectangular grid database
  virtual
  void Add(const std::vector<double>& data, const std::string& name) = 0;

  //! Write output file
  virtual
  void Write(const std::string& name) = 0;
};

typedef std::shared_ptr<ThreeDGrid> ThreeDGridPtr;

#endif
