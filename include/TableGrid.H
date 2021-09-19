#ifndef _TableGrid_H
#define _TableGrid_H

#include <sstream>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/multi_array.hpp>

#include <ThreeDGrid.H>

class TableGrid : public ThreeDGrid
{
protected:
  // Knots
  std::vector<float> XX, YY, ZZ;

  // Data
  typedef boost::multi_array<float, 3> dtype;
  typedef boost::shared_ptr<dtype> dptr;
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

typedef boost::shared_ptr<TableGrid> TableGridPtr;

#endif
