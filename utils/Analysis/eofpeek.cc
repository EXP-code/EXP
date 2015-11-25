#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

int main (int nc, char** nv)
{
  if (nc != 2) {
    std::cout << "Usage: " << nv[0] << " <eof cache>" << std::endl;
  }

  std::ifstream in(nv[1]);
  if (!in) {
    std::cerr << "EmpCylSL::cache_grid: error opening file named <" 
	      << nv[1] << ">" << std::endl;
    return 1;
  }


  int MMAX, NUMX, NUMY, NMAX, NORDER, tmp;
  double RMIN, RMAX, ASCALE, HSCALE;
  bool DENS, CMAP;

  in.read((char *)&MMAX,   sizeof(int));
  in.read((char *)&NUMX,   sizeof(int));
  in.read((char *)&NUMY,   sizeof(int));
  in.read((char *)&NMAX,   sizeof(int));
  in.read((char *)&NORDER, sizeof(int));
  in.read((char *)&tmp,    sizeof(int)); 
  if (tmp) DENS = true; else DENS = false;
  in.read((char *)&tmp,    sizeof(int)); 
  if (tmp) CMAP = true; else CMAP = false;
  in.read((char *)&RMIN,   sizeof(double));
  in.read((char *)&RMAX,   sizeof(double));
  in.read((char *)&ASCALE, sizeof(double));
  in.read((char *)&HSCALE, sizeof(double));
  
  // BEGIN header
  std::ostringstream sout;
  sout << "----- Cylindrical parameters read from <" << nv[1] << "> ";

  std::cout << std::setfill('-') << std::setw(70) << '-' << std::endl
	    << std::setw(70) << std::left << sout.str()  << std::endl
	    << std::setw(70) << '-' << std::setfill(' ') << std::endl;
  // END header

  std::cout << std::setw(12) << std::left << "MMAX"   << MMAX   << std::endl;
  std::cout << std::setw(12) << std::left << "NUMX"   << NUMX   << std::endl;
  std::cout << std::setw(12) << std::left << "NUMY"   << NUMY   << std::endl;
  std::cout << std::setw(12) << std::left << "NMAX"   << NMAX   << std::endl;
  std::cout << std::setw(12) << std::left << "NORDER" << NORDER << std::endl;
  std::cout << std::setw(12) << std::left << "CMAP"   << (CMAP ? "true" : "false") << std::endl;
  std::cout << std::setw(12) << std::left << "DENS"   << (DENS ? "true" : "false") << std::endl;
  std::cout << std::setw(12) << std::left << "RMIN"   << RMIN   << std::endl;
  std::cout << std::setw(12) << std::left << "RMAX"   << RMAX   << std::endl;
  std::cout << std::setw(12) << std::left << "ASCALE" << ASCALE << std::endl;
  std::cout << std::setw(12) << std::left << "HSCALE" << HSCALE << std::endl;

  std::cout << std::setfill('-') << std::left << std::setw(70) << '-' << std::endl;

  return 0;
}
