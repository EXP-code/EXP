#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <memory>
#include <string>

#include <Eigen/Eigen>

#include <localmpi.H>

// Structure from SLGridSph
//
class TableSph 
{
 public:
  int l;

  Eigen::VectorXd ev;
  Eigen::MatrixXd ef;
};


int main(int argc, char **argv)
{
  if (argc<2) {
    std::cerr << "Usage: " << argv[0] << " <cache file>" << std::endl;
    exit(-1);
  }

  
  std::ifstream in(argv[1]);
  if (!in) {
    std::cerr << "Could not read cache file <" << argv[1] << ">"
	      << std::endl;
  }

  int lmax, nmax, numr, cmap;
  double rmin, rmax, scl;

  if (myid==0) 
    std::cerr << "SLGridSph::read_cached_table: trying to read cached table . . ."
	      << std::endl;

  in.read((char *)&lmax, sizeof(int));
  in.read((char *)&nmax, sizeof(int));
  in.read((char *)&numr, sizeof(int));
  in.read((char *)&cmap, sizeof(int));
  in.read((char *)&rmin, sizeof(double));
  in.read((char *)&rmax, sizeof(double));
  in.read((char *)&scl,  sizeof(double));

  std::cout << std::left
	    << std::setw(20) << "# Lmax"
	    << std::setw(20) << nmax << std::endl
	    << std::setw(20) << "# Nmax"
	    << std::setw(20) << lmax << std::endl
	    << std::setw(20) << "# NumR"
	    << std::setw(20) << numr << std::endl
	    << std::setw(20) << "# cmap"
	    << std::setw(20) << cmap << std::endl
	    << std::setw(20) << "# r_min"
	    << std::setw(20) << rmin << std::endl
	    << std::setw(20) << "# r_max"
	    << std::setw(20) << rmax << std::endl
	    << std::setw(20) << "# r_scl"
	    << std::setw(20) << scl << std::endl
	    << "#" << std::endl << std::right;
    

  std::vector<TableSph> table(lmax+1);

  for (int l=0; l<=lmax; l++) {

    in.read((char *)&table[l].l, sizeof(int));

				// Double check
    if (table[l].l != l) {
      if (myid==0)
	std::cerr << "SLGridSph: error reading <" << argv[1] << ">" << std::endl
		  << "SLGridSph: l: read value (" << table[l].l 
		  << ") != internal value (" << l << ")" << std::endl;
	return 0;
    }

    table[l].ev.resize(nmax);
    table[l].ef.resize(nmax, numr);

    for (int j=0; j<nmax; j++)
      for (int i=0; i<numr; i++)
	in.read((char *)&table[l].ef(j, i), sizeof(double));

    Eigen::MatrixXd tmp(numr, nmax);
    in.read((char *)table[l].ev.data(), nmax*sizeof(double));
    in.read((char *)tmp.data(), numr*nmax*sizeof(double));
    table[l].ef = tmp.transpose();
  }

  std::cerr << "Success reading cache file!!" << std::endl;

  for (int l=0; l<=lmax; l++) {
    std::cout << "# l=" << l << std::endl;
    for (int i=0; i<numr; i++) {
      std::cout << std::setw( 5) << i;
      for (int j=0; j<=std::min(nmax, 4); j++)
	std::cout << std::setw(20) << table[l].ef(j, i);
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  return 0;
}
