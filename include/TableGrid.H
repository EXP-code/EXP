#ifndef _TableGrid_H
#define _TableGrid_H

#include <sstream>
#include <vector>
#include <memory>
#include <map>

#include <memory>

#include "ThreeDGrid.H"

/**
   This implementation of ThreeDGrid will write ascii files with
   records of the form X, Y, Z, D1, D2, D3, ... where Dn are the data
   fields loaded by n called to the Add() member.
 */
class TableGrid : public ThreeDGrid
{
protected:
  // Knots
  std::vector<float> XX, YY, ZZ;

  // Data
  typedef Dynamic3dArray<double> dtype;
  typedef std::shared_ptr<dtype> dptr;
  std::map<std::string, dptr> arrays;

public:
  //! Constructor
  TableGrid(int nx, int ny, int nz,
	    double xmin, double xmax,
	    double ymin, double ymax,
	    double zmin, double zmax);

  //! Add data
  void Add(const std::vector<double>& data, const std::string& name);

  //! Write output file
  void Write(const std::string& name);
};

typedef std::shared_ptr<TableGrid> TableGridPtr;

#endif
