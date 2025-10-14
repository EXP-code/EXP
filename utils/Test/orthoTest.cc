/*****************************************************************************
 *  Description:
 *  -----------
 *
 *  Computes orthogonal function expansion from the biorthogonal basis set
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
 *  revision MDW 09/04/2025
 *
 ***************************************************************************/

#include <filesystem>
#include <functional>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <string>
#include <cstring>
#include <vector>
#include <stdexcept>		// For potential error handling
#include <type_traits>		// For std::enable_if and type traits

#include <Eigen/Eigen>

#include <localmpi.H>
#include <Progress.H>
#include <cxxopts.H>
#include <Timer.H>

#include <massmodel.H>
#include <isothermal.H>
#include <model3d.H>
#include <biorth.H>
#include <SLGridMP2.H>
#include <Biorth2Ortho.H>

#ifdef FPETRAP
#include <fenv.h>
#endif

// Define floating point type and matrix/vector types
//
using Real = long double;
using MatrixXld = Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>;
using VectorXld = Eigen::Vector<Real, Eigen::Dynamic>;

// Overload output operator for MatrixXld, needed for __float128
//
std::ostream& operator<<(std::ostream& os, const MatrixXld& M)
{
  for (int i=0; i<M.rows(); i++) {
    for (int j=0; j<M.cols(); j++) {
      os << std::setw(18) << static_cast<double>(M(i, j));
    }
    os << std::endl;
  }
  return os;
};


int
main(int argc, char **argv)
{
  // MPI initialization
  //
  local_init_mpi(argc, argv);

  // Parameters
  //
  int L, M, lmax, nmax, MODEL, NUMDF, NRGRID, DIVERGE;
  int type, NUMR, ng;
  double DIVERGE_RFAC, rmin, rmax, RA, scl;
  std::string INFILE, PREFIX, suffix;
  BiorthFcts3d TYPE;
  
  // Parse command line
  //
  std::string message = "\nComputes biorthogonal to orthogonal transformation\n";

  cxxopts::Options options(argv[0], message);

  options.add_options()
    ("h,help", "Print this help message")
    ("potential", "Use potential expansion (default is density)")
    ("alternate", "Use both density and potential expansion to contruct new basis (default is density only)")
    ("sum", "Combine density and potential to contruct new basis (default is density")
    ("laguerre", "Use Laguerre polynomials for Gram-Schmidt orthogonalization (default is to use biorthogonal functions)")
    ("classic", "Use classical Gram-Schmidt orthogonalization (default is modified Gram-Schmidt)")
    ("nobar", "Do not show progress bar")
    ("weight", "Use density weighting for orthonormal functions")
    ("p,progress", "Show progress bar")
    ("d,diag", "Print diagnostic output")
    ("l,ll", "Harmonic order (only used for three-dimensional spherical models)",
     cxxopts::value<int>(L)->default_value("1"))
    ("m,mm", "Harmonic order (only used for two-dimensional cylidrical models)",
     cxxopts::value<int>(M)->default_value("1"))
    ("L,Lmax", "Maximum harmonic order",
     cxxopts::value<int>(lmax)->default_value("4"))
    ("n,Nmax", "Maximum radial order",
     cxxopts::value<int>(nmax)->default_value("12"))
    ("T,TYPE", "Biothogonal basis type (bessel[0], clutton_brock[1], hernquist[2], sturm[3])",
     cxxopts::value<int>(type)->default_value("3"))
    ("MODEL", "Model type",
     cxxopts::value<int>(MODEL)->default_value("0"))
    ("NUMDF", "number of distribution function grid points",
     cxxopts::value<int>(NUMDF)->default_value("400"))
    ("rmin", "minimum expansion radius",
     cxxopts::value<double>(rmin)->default_value("-1.0"))
    ("rmax", "maximum expansion radius",
     cxxopts::value<double>(rmax)->default_value("-1.0"))
    ("NUMR", "Number of radial grid points for SLGridSph",
     cxxopts::value<int>(NUMR)->default_value("200"))
    ("scl", "Coordinate mapping scale for SLGridSph",
     cxxopts::value<double>(scl)->default_value("0.05"))
    ("a,RA", "Osipkov-Merritt radius",
     cxxopts::value<double>(RA)->default_value("1.0e20"))
    ("DIVERGE_RFAC", "Inner power-law slope",
     cxxopts::value<double>(DIVERGE_RFAC)->default_value("0"))
    ("DIVERGE", "Use inner power-law slope extrapolation",
     cxxopts::value<int>(DIVERGE)->default_value("0"))
    ("NRGRID", "Grid points for BiorthGrid",
     cxxopts::value<int>(NRGRID)->default_value("512"))
    ("NG", "Grid points for scalar products",
     cxxopts::value<int>(ng)->default_value("512"))
    ("i,INFILE", "Model file",
     cxxopts::value<std::string>(INFILE)->default_value("model.dat"))
    ("o,PREFIX", "Output file prefix",
     cxxopts::value<std::string>(PREFIX)->default_value("orthoTest"))
    ;

  cxxopts::ParseResult vm;

  try {
    vm = options.parse(argc, argv);
  } catch (cxxopts::OptionException& e) {
    if (myid==0) std::cout << "Option error: " << e.what() << std::endl;
    MPI_Finalize();
    return -1;
  }
  
  // Get help
  //
  if (vm.count("help")) {
    if (myid==0) std::cout << std::endl << options.help() << std::endl;
    MPI_Finalize();
    return 0;
  }

  // Assign biorthogonal type
  //
  switch(type) {
  case 0:
    TYPE = bessel;
    break;
  case 1:
    TYPE = clutton_brock;
    break;
  case 2:
    TYPE = hernquist;
    break;
  case 3:
    TYPE = sturm;
    break;
  default:
    std::cerr << "No such biorthonal TYPE: " << type << std::endl;
    return -1;
  }

#ifdef FPETRAP
  fpeinit(0);
#endif


  // Initilize model and generate Osipkov-Merritt distribution function
  //
  SphModTblPtr m;
  AxiSymModPtr model;

  switch (MODEL) {
  case file:
    m = std::make_shared<SphericalModelTable>(INFILE, DIVERGE, DIVERGE_RFAC);
    m->setup_df(NUMDF, RA);
    model = m;
    Model3dNames[0] = INFILE;	// Assign filename to ID string
    break;

  case isothermal:
    model = std::make_shared<IsothermalSphere>(1.0, 1.0e-3, 1.0e3);
    break;

  default:
    std::cerr << "Illegal model: " << MODEL << std::endl;
    exit(-1);
  }

  if (rmin < 0.0) rmin = model->get_min_radius();
  if (rmax < 0.0) rmax = model->get_max_radius();

  
  // Initilize biorthogonal functions 
  //
  std::string cache= ".slgrid_sph_cache";

  // SLGridSph::cache = 0;		// Turn off cache reading

  AxiSymBioPtr biorth;
  switch (TYPE) {
  case bessel:
    biorth = std::make_shared<BSSphere>(rmax, nmax, lmax);
    break;
  case clutton_brock:
    biorth = std::make_shared<CBSphere>();
    break;
  case hernquist:
    biorth = std::make_shared<HQSphere>();
    break;
  case sturm:
    biorth = std::make_shared<Sturm>(m, lmax, nmax, NUMR, cache, scl);
    break;
  default:
    std::cerr << "Illegal type: " << TYPE << std::endl;
    exit(-1);
  }

  bool weight = vm.count("weight")>0;

  Biorth2Ortho trans(biorth, lmax, nmax, NRGRID, rmin, rmax, scl, weight);

  if (vm.count("nobar")) trans.noBar();
  trans.useLaguerre(vm.count("laguerre")>0);
  trans.useAlternate(vm.count("alternate")>0);
  trans.useSum(vm.count("sum")>0);
  trans.useClassic(vm.count("classic")>0);
  trans.useDensity(vm.count("potential")==0);
  
  trans.generate();
  trans.output(PREFIX);
  trans.writeH5(PREFIX + ".h5");

  // All done!
  //
  MPI_Finalize();

  return 0;
}
