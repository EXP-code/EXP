// Helper routine to check an orthogonality matrix and throw an
// exception if tolerance is not met

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include <Eigen/Eigen>		// For MatrixXd

#include <omp.h>		// For OpenMP API

#include <localmpi.H>		// MPI support
#include <libvars.H>		// For orthoTol

std::tuple<bool, double, std::vector<double>>
orthoCompute(const std::vector<Eigen::MatrixXd>& tests)
{
  // Number of possible threads
  int nthrds = omp_get_max_threads();

  // Worst so far
  std::vector<double> worst(nthrds), lworst(tests.size());;

  // Rank
  int nmax = tests[0].rows();
  
  // Test loop
  for (int l=0; l<tests.size(); l++) {

    // Initialize test array
    std::fill(worst.begin(), worst.end(), 0.0);

#pragma omp parallel for
    for (int nn=0; nn<nmax*nmax; nn++) {
      int tid = omp_get_thread_num();
      int  n1 = nn/nmax;
      int  n2 = nn - n1*nmax;
      
      if (n1==n2)
	worst[tid] = std::max<double>(worst[tid],
				      fabs(1.0 - tests[l](n1, n2)));
      else
	worst[tid] = std::max<double>(worst[tid],
				      fabs(tests[l](n1, n2)));
    }
    // END: unrolled loop
    
    lworst[l] = *std::max_element(worst.begin(), worst.end());
  }
  // END: harmonic order loop

  double worst_ever = *std::max_element(lworst.begin(), lworst.end());

  if (worst_ever > __EXP__::orthoTol) {
    return {false, worst_ever, lworst};
  } else {
    return {true, worst_ever, lworst};
  }
}


void orthoCompare(const std::vector<Eigen::MatrixXd>& tests,
		  const std::string& classname, const std::string& indexname)
{
  auto [good, worst, lworst] = orthoCompute(tests);

  if (good) {
    if (myid==0) 
      std::cout << classname + ": biorthogonal check passed" << std::endl;
  } else {
    std::cout << classname << ": orthogonality failure" << std::endl
	      << std::right
	      << std::setw(4) << indexname
	      << std::setw(16) << "Worst" << std::endl;
    for (int l=0; l<lworst.size(); l++) {
      std::cout << std::setw(4) << l << std::setw(16) << lworst[l] << std::endl;
    }
  }
}
