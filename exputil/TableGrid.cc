#include <iomanip>
#include <fstream>

#include <TableGrid.H>

TableGrid::TableGrid(int nx, int ny, int nz,
		 double xmin, double xmax,
		 double ymin, double ymax,
		 double zmin, double zmax) :
  ThreeDGrid(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax)
{
  // Set knots
  //
  XX.resize(nx);
  YY.resize(ny);

  if (nz<0) nz = 1;
  if (nz>1) ZZ.resize(nz);

  for (int i=0; i<nx; i++) {
    XX[i] = xmin + (xmax - xmin)*i/(nx-1);
  }

  for (int i=0; i<ny; i++) {
    YY[i] = ymin + (ymax - ymin)*i/(ny-1);
  }

  if (nz>1) {
    for (int i=0; i<nz; i++)
      ZZ[i] = zmin + (zmax - zmin)*i/(nz-1);
  }    

}

void TableGrid::Add(const std::vector<double>& data, const std::string& name)
{
  // Set array size
  //
  auto D = std::make_shared<dtype>(nx, ny, nz);

  // Insert grid data
  //
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      for (int k=0; k<nz; k++) {
	(*D)[i][j][k] = data[(i*ny+j)*nz + k];
      }
    }
  }

  // Add to database
  //
  arrays[name] = D;
}


void TableGrid::Write(const std::string& name)
{
  // Create the filename with the correct extension for Paraview
  //
  std::ostringstream filename;
  filename << name << ".ascii";

  std::ofstream out(filename.str());
  
  // Make header
  //
  out << "# [ 1] XX" << std::endl
      << "# [ 2] YY" << std::endl
      << "# [ 3] ZZ" << std::endl;
  int cnt = 4;
  for (auto & p : arrays) {
    out << "# [" << std::setw(2) << cnt++ << "] " << p.first << std::endl;
  }
  out << std::endl;


  // Insert grid data
  //
  for (int i=0; i<nx; i++) {
    for (int j=0; j<ny; j++) {
      for (int k=0; k<nz; k++) {
	out << std::setw(18) << XX[i]
	    << std::setw(18) << YY[j]
	    << std::setw(18) << ZZ[k];
	for (auto & p : arrays)
	  out << std::setw(18) << (*p.second)[i][j][k];
	out << std::endl;
      }
      if (nz>1) out << std::endl;
    }
    out << std::endl;
  }
}
