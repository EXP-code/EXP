/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Read in coefficients and compute gnuplot slices, and compute
 *  volume for rendering
 *
 *
 *  Call sequence:
 *  -------------
 *
 *  Parameters:
 *  ----------
 *
 *
 *  Returns:
 *  -------
 *
 *
 *  Notes:
 *  -----
 *
 *
 *  By:
 *  --
 *
 *  MDW 11/28/08
 *
 ***************************************************************************/

				// C++/STL headers
#include <filesystem>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>

using namespace std;

                                // System libs
#include <sys/time.h>
#include <sys/resource.h>

				// MDW classes
#include <numerical.H>
#include <ParticleReader.H>
#include <interp.H>
#include <massmodel.H>
#include <SphereSL.H>
#include <DataGrid.H>

// Library support
#include <localmpi.H>
#include <foarray.H>
#include <EXPini.H>
//#include <global.H>
#include <KDtree.H>
//#include <libvars.H>

//using namespace __EXP__;

// generic header includes: requires both pybind11 and eigen
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

// eigen includes
#include <Eigen/StdVector>
#include <Eigen/Dense>
using Eigen::MatrixXd;

// Globals
//
std::string OUTFILE;
double RMIN, RMAX;
int OUTR, LMAX, NMAX, L1, L2, N1, N2;
bool VOLUME, SURFACE, PROBE;

// Center offset
//
std::vector<double> c0 = {0.0, 0.0, 0.0};



Eigen::MatrixXd
makesl(std::vector<double> mass,
       std::vector<double> xpos, std::vector<double> ypos, std::vector<double> zpos,
       string MODFILE,
       int LMAX=4,
       int NMAX=10)
{

  // count the number of particles
  int nparticles = xpos.size();

  auto halo = std::make_shared<SphericalModelTable>(MODFILE);

  SphereSL::mpi  = false;
  SphereSL::NUMR = 4000;
  //SphereSL::HEXP = 0.002;//Hexp;
  double rscale  = 1.;

  SphereSL ortho(halo, LMAX, NMAX, 1, rscale, true);

  // set up for accumulation
  ortho.reset_coefs();

  // loop through particles
  for (int i=0;i<nparticles;i++) {
        ortho.accumulate(xpos[i],ypos[i],zpos[i],mass[i]);
      }

  std::cout << "done" << std::endl;

      //------------------------------------------------------------

  std::cout << "Making coefficients . . . " << std::flush;
  ortho.make_coefs();
  std::cout << "done" << std::endl;

  // reshape the ortho.expcoef into a MatrixXd
  // allocate workspace
  MatrixXd outcoefs;
  outcoefs = ortho.retrieve_coefs();

  return outcoefs;
}


PYBIND11_MODULE(simpleSL, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("makesl", &makesl, "Make an SL");
}
