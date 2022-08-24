#include <algorithm>
#include <iomanip>
#include <sstream>

#include <FieldGenerator.H>

FieldGenerator::FieldGenerator(const std::vector<double> &time,
			       const std::vector<double> &pmin,
			       const std::vector<double> &pmax,
			       const std::vector<int>    &grid
			       ) :
  time(time), pmin(pmin), pmax(pmax), grid(grid)
{
  // Nothing so far
}

void FieldGenerator::check_times(CoefPtr coefs)
{
  std::vector<double> ctimes = coefs->Times();
  std::sort(ctimes.begin(), ctimes.end());

  for (auto t : time) {
    if (std::find(ctimes.begin(), ctimes.end(), t) == ctimes.end()) {
      std::ostringstream sout;
      sout << "FieldGenerator: requested time <" << t << "> "
	   << "not in DB" << std::endl;
      throw std::runtime_error(sout.str());
   }
  }
}

std::map<double, std::map<std::string, Eigen::MatrixXd>>
FieldGenerator::slices(BasisPtr basis, CoefPtr coefs)
{
  std::map<double, std::map<std::string, Eigen::MatrixXd>> ret;
  return ret;
}

void FieldGenerator::vtk_slices(BasisPtr basis, CoefPtr coefs)
{
}


std::map<double, std::map<std::string, Eigen::Tensor<float, 3>>>
FieldGenerator::volumes(BasisPtr basis, CoefPtr coefs)
{
  std::map<double, std::map<std::string, Eigen::Tensor<float, 3>>> ret;
  return ret;
}

void FieldGenerator::vtk_volumes(BasisPtr basis, CoefPtr coefs)
{
}


