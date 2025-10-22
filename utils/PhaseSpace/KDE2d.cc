#include <KDE2d.H>

namespace KDE
{
  // Function to compute the Gaussian kernel in pixel units
  //
  std::vector<double> KDE2d::gaussian_kernel_2d(int xsize, int ysize)
  {
    std::vector<double> kernel(xsize * ysize);
    double sum = 0.0;
    for (int i = -xsize / 2; i < xsize / 2 + 1; i++) {

      double dx = delx*i/numx;
      double facx = exp(-0.5*dx*dx/(sigmax*sigmax));
  
      for (int j = -ysize / 2; j < ysize / 2 + 1; j++) {

	double dy = dely*j/numy;
	double facy = exp(-0.5*dy*dy/(sigmay*sigmay)); 

	kernel[(i + xsize / 2) * ysize + (j + ysize / 2)] = facx * facy;
	
	sum += facx * facy;
      }
    }
    
    // Normalize the kernel
    //
    for (auto & v : kernel) v /= sum;

    return kernel;
  }

  
  void KDE2d::kde_fft_2d()
  {
    // Sanity check
    //
    if (numx != grid.rows() || numy != grid.cols())
      throw std::runtime_error("Invalid grid dimensions");

    // Compute the Gaussian kernel
    //
    int kern_xsize = std::max<int>(minKsize,
				   std::floor(sigmax/delx*numx*minKsize));

    if (kern_xsize % 2 == 0) kern_xsize += 1; // Make the size odd
    
    int kern_ysize = std::max<int>(minKsize,
				   std::floor(sigmay/dely*numy*minKsize));

    if (kern_ysize % 2 == 0) kern_ysize += 1; // Make the size odd
    
    std::vector<double> kern = gaussian_kernel_2d(kern_xsize, kern_ysize);
    
    if (debug) {
      std::ofstream fout("kern0.dat");
      fout << kern_xsize << " " << kern_ysize << std::endl;
      for (int j=0; j<kern_ysize; j++) {
	for (int i=0; i<kern_xsize; i++) {
	  fout << std::setw( 8) << i 
	       << std::setw( 8) << j
	       << std::setw(18) << kern[i*kern_ysize + j] << std::endl;
	}
      }
    }
    
    // Pad the grid and kernel to the same size
    //
    int padded_xsize = numx + kern_xsize - 1;
    int padded_ysize = numy + kern_ysize - 1;
    
    std::vector<double> padded_grid(padded_xsize * padded_ysize, 0.0);
    std::vector<double> padded_kern(padded_xsize * padded_ysize, 0.0);
    
    for (int i=0; i<numx; i++) {
      for (int j=0; j<numy; j++) {
	padded_grid[i*padded_ysize + j] = grid(i, j);
      }
    }
    

    for (int i=0; i<kern_xsize; i++) {

      int ii = i - kern_xsize/2;
      if (ii < 0) ii += padded_xsize;
      
      for (int j=0; j<kern_ysize; j++) {
	
	int jj = j - kern_ysize/2;
	if (jj < 0) jj += padded_ysize;
	
	padded_kern[ii*padded_ysize + jj] = kern[i*kern_ysize + j];
      }
    }

    if (debug) {
      std::ofstream kdat("kern.dat");
      std::ofstream ddat("data.dat");
      std::ofstream gdat("grid.dat");
      kdat << padded_xsize << " " << padded_ysize << std::endl;
      ddat << padded_xsize << " " << padded_ysize << std::endl;
      gdat << numx << " " << numy << std::endl;
      for (int j=0; j<padded_ysize; j++) {
	for (int i=0; i<padded_xsize; i++) {
	  kdat << std::setw(18) << delx*i/numx
	       << std::setw(18) << dely*j/numy
	       << std::setw(18) << padded_kern[i*padded_ysize + j] << std::endl;
	  ddat << std::setw(18) << xmin + delx*i/numx
	       << std::setw(18) << ymin + dely*j/numy
	       << std::setw(18) << padded_grid[i*padded_ysize + j] << std::endl;
	}
      }
      for (int j=0; j<numy; j++) {
	for (int i=0; i<numx; i++) {
	  gdat << std::setw(18) << xmin + delx*i/numx
	       << std::setw(18) << ymin + dely*j/numy
	       << std::setw(18) << grid(i, j) << std::endl;
	}
      }
    }

    // Perform FFT on the padded grid and kernel using fftw C arrays
    //
    int asize = padded_xsize * padded_ysize;
    size_t sz = sizeof(fftw_complex) * asize;
    
    auto grid_fft = (fftw_complex*)fftw_malloc(sz);
    auto kern_fft = (fftw_complex*)fftw_malloc(sz);
    auto rslt_fft = (fftw_complex*)fftw_malloc(sz);
    
    for (int i = 0; i < asize; i++) {
      grid_fft[i][0] = padded_grid[i];
      grid_fft[i][1] = 0.0;
      kern_fft[i][0] = padded_kern[i];
      kern_fft[i][1] = 0.0;
    }
  
    auto plan_forward_grid =
      fftw_plan_dft_2d(padded_xsize, padded_ysize, grid_fft, grid_fft,
		       FFTW_FORWARD, FFTW_ESTIMATE);

    auto plan_forward_kern =
      fftw_plan_dft_2d(padded_xsize, padded_ysize, kern_fft, kern_fft,
		       FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_execute(plan_forward_grid);
    
    fftw_execute(plan_forward_kern);

    // Multiply the FFTs to compute the convolution
    //
    for (int i=0; i<asize; i++) {
      rslt_fft[i][0] = grid_fft[i][0] * kern_fft[i][0] - grid_fft[i][1] * kern_fft[i][1];
      rslt_fft[i][1] = grid_fft[i][0] * kern_fft[i][1] + grid_fft[i][1] * kern_fft[i][0];
    }
    
    // Perform inverse FFT
    //
    auto plan_backward =
      fftw_plan_dft_2d(padded_xsize, padded_ysize, rslt_fft, rslt_fft,
		       FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute(plan_backward);

    // Copy the result to remove the padding
    //
    smooth.resize(numx, numy);
    for (int i=0; i<numx; i++) {
      for (int j=0; j<numy; j++) {
	smooth(i, j) = rslt_fft[i*padded_ysize + j][0] / asize;
      }
    }
    
    if (debug) {
      std::ofstream rdat("rslt.dat");
      std::ofstream odat("outp.dat");
      rdat << padded_xsize << " " << padded_ysize << std::endl;
      odat << numx << " " << numy << std::endl;
      for (int j=0; j<padded_ysize; j++) {
	for (int i=0; i<padded_xsize; i++) {
	  rdat << std::setw(18) << xmin + delx*i/numx
	       << std::setw(18) << ymin + dely*j/numy
	       << std::setw(18) << rslt_fft[i*padded_ysize + j][0] << std::endl;
	}
      }
      for (int j=0; j<numy; j++) {
	for (int i=0; i<numx; i++) {
	  odat << std::setw(18) << xmin + delx*i/numx
	       << std::setw(18) << ymin + dely*j/numy
	       << std::setw(18) << smooth(i, j) << std::endl;
	}
      }
    }
    
    // Clean up
    //
    fftw_destroy_plan(plan_forward_grid);
    fftw_destroy_plan(plan_forward_kern);
    fftw_destroy_plan(plan_backward);
  
    fftw_free(grid_fft);
    fftw_free(kern_fft);
    fftw_free(rslt_fft);
  }

  // Function to perform 2D KDE using FFT
  //
  void
  KDE2d::grid_pairs(const std::vector<std::pair<double, double>>& data)
  {
    // Sanity check
    //
    if (data.size() == 0)
      throw std::invalid_argument("data must be non-empty");

    // Create a grid to evaluate KDE
    //
    grid.resize(numx, numy);
    grid.setZero();
    for (const auto& point : data) {
      int x = static_cast<int>( (point.first  - xmin)/delx * numx);
      int y = static_cast<int>( (point.second - ymin)/dely * numy);

      if (x >= 0 && x < numx && y >= 0 && y < numy) {
	grid(x, y) += 1.0;
      }
      else {
	std::cout << "Out of range: " << point.first << " " << point.second << std::endl;
      }
    }
  }

  void
  KDE2d::grid_array(const std::vector<double>& X, const std::vector<double>& Y)
  {
    // Sanity check
    //
    if (X.size() != Y.size())
      throw std::invalid_argument("x and y must be the same size");

    if (X.size() == 0)
      throw std::invalid_argument("data must be non-empty");

    // Create a grid to evaluate KDE
    //
    grid.resize(numx, numy);
    grid.setZero();
    for (int i=0; i<X.size(); i++) {
      int x = std::floor( (X[i] - xmin)/delx * numx);
      int y = std::floor( (Y[i] - ymin)/dely * numy);
      if (x >= 0 && x < numx && y >= 0 && y < numy) {
	grid(x, y) += 1.0;
      }
      else {
	std::cout << "Out of range: " << X[i] << " " << Y[i] << std::endl;
      }
    }
  }

}
// END namespace KDE
