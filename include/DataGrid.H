#ifndef _DataGrid_H
#define _DataGrid_H

#include <vector>
#include <sstream>
#include <memory>

#include <config_exp.h>

#include <TableGrid.H>
#include <VtkGrid.H>

/**
   This implementation of ThreeDGrid instantiates at VtkGrid if VTK is
   available and TableGrid if not.

   This is really syntactic glue that allows one call without checking
   for VTK in application codes.
*/
class DataGrid : public ThreeDGrid
{
private:

  ThreeDGridPtr ptr;

public:

  //! Constructor
  DataGrid(int nx, int ny, int nz,
	   double xmin, double xmax,
	   double ymin, double ymax,
	   double zmin, double zmax,
	   std::string type="VTK")
  {
    std::transform(type.begin(), type.end(), type.begin(), 
                   [](unsigned char c){ return std::tolower(c); });

    if (type.compare("vtk") == 0)
      ptr = std::make_shared<VtkGrid>(nx, ny, nz,
				      xmin, xmax,
				      ymin, ymax,
				      zmin, zmax);
    else
      ptr = std::make_shared<TableGrid>(nx, ny, nz,
					xmin, xmax,
					ymin, ymax,
					zmin, zmax);
  }
    
  //! Add data
  void Add(const std::vector<double>& data, const std::string& name)
  {
    ptr->Add(data, name);
  }

  //! Write output file
  void Write(const std::string& name)
  {
    ptr->Write(name);
  }
};

typedef std::shared_ptr<DataGrid> DataGridPtr;

#endif
