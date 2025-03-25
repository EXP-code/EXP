#ifndef KDE2D_H
#define KDE2D_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include <Eigen/Eigen>

#include <fftw3.h>

#include <cxxopts.H>

namespace KDE
{
  class KDE2d
  {
  protected:

    //@{
    //! Parameters
    bool debug = false;
    const int minKsize = 7;
    int numx, numy;
    double xmin, xmax, ymin, ymax, sigmax, sigmay, delx, dely;
    //@}

    //! Create a smooth 2D gaussian kernel
    std::vector<double>
    gaussian_kernel_2d(int xsize, int ysize);
  
    //! Do the FFT convolution
    void kde_fft_2d();

    //@{
    //! Make 2D grid
    void grid_pairs(const std::vector<std::pair<double, double>>& data);
    void grid_array(const std::vector<double>& x, const std::vector<double>& y);
    //@}

    //! Eigen matrices
    Eigen::MatrixXd grid, smooth;

  public:

    //! Constructor
    KDE2d(int numx, int numy,
	  double xmin, double xmax,
	  double ymin, double ymax,
	  double sigmax, double sigmay) :
      numx(numx), numy(numy), xmin(xmin), xmax(xmax),
      ymin(ymin), ymax(ymax), sigmax(sigmax), sigmay(sigmay)
    {
      delx = xmax - xmin;
      dely = ymax - ymin;
      if (delx<=0 || dely<=0)
	throw std::invalid_argument("Grid size must be positive");
    }

    //! Get the smoothed density from an input grid
    const Eigen::MatrixXd& operator()(const Eigen::MatrixXd& GRID)
    {
      if (GRID.rows() != numx || GRID.cols() != numy)
	throw std::invalid_argument("Grid size must match KDE2d object");
      grid = GRID;
      kde_fft_2d();
      return smooth;
    }

    //! Get the smoothed density from an input data pairs
    const Eigen::MatrixXd& operator()(const std::vector<std::pair<double, double>>& data)
    {
      grid_pairs(data);
      kde_fft_2d();
      return smooth;
    }

    //! Get the smoothed density from an input data arrays
    const Eigen::MatrixXd& operator()
      (const std::vector<double>& x, const std::vector<double>& y)
    {
      grid_array(x, y);
      kde_fft_2d();
      return smooth;
    }

    //! Set debug mode
    void setDebug() { debug = true; }
  };
}

#endif // KDE2D_H
