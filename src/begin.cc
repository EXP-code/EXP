/*
  Get the whole thing going by:
  -----------------------------
  1) reading in initial phase space
  2) computing potential initially
  3) dumping first phase-space and log file
  4) write parameter file
*/

#include <sys/types.h>
#include <unistd.h>		// For getpid

#include "expand.H"
#include "ExternalCollection.H"
#include "OutputContainer.H"

void begin_run(void)
{
  //===================================
  // Initialize cuda device(s)
  //===================================

  initialize_cuda();

  //===================================
  // Make the instance containers
  //===================================
  
  comp     = new ComponentContainer;
  external = new ExternalCollection;
  output   = new OutputContainer;

  //===================================
  // Make the Vector class mutex
  //===================================
  
  int errcode = pthread_mutex_init(&mem_lock, NULL), nOK = 0;
  if (errcode) {
    cerr << "Process " << myid 
	 << ": failed to make the Vector class memory lock, errcode=" 
	 << errcode << endl;
    nOK = 1;
  }

  MPI_Allreduce(&nOK, &errcode, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if (errcode) {
    MPI_Finalize();
    exit(115);
  }

  //===================================
  // Make the kick/drift thread vectors
  //===================================
  
  posvel_data = vector<thrd_pass_posvel>(nthrds);
  posvel_thrd = vector<pthread_t>       (nthrds);

  //==============================
  // Initialize multistepping
  //==============================

  initialize_multistep();

  //===================================
  // Initialize phase-space components 
  //===================================

  comp->initialize();

  //===================================
  // Initialize external forces, if any
  //===================================

  external->initialize();
  
  //===============================
  // Compute initial acceleration
  //===============================

  initializing = true;
  
  //================================
  // Multistep level initialization
  //================================

  if (multistep) {

    comp->multistep_reset();

    //=====================================
    // Compute coefficients at every level
    //=====================================

    for (int M=0; M<=multistep; M++) comp->compute_expansion(M);
    //                          ^
    //                          |
    // Loop on all levels-------+

    //========================
    // Compute full potential
    //========================

    comp->compute_potential(0);
    //                      ^
    //                      |
    //   All time levels----/

    //==============================
    // Initialize multistep levels
    //==============================

    adjust_multistep_level();
    //  ^
    //  |  
    //  \----All on first call
    
  }

  //===========================================
  // Compute coefficients (again if multistep)
  //===========================================

  if (multistep) comp->multistep_reset();

  for (int M=0; M<=multistep; M++) comp->compute_expansion(M);

  comp->compute_potential(0);

  initializing = false;

  //===================================
  // Initialize output routines
  //===================================

  output->initialize();
  output->Run(0);

  //======================
  // Write parameter file 
  //======================

  write_parm();

}


void initialize_cuda()
{
#if HAVE_LIBCUDA==1
  int deviceCount = 0;

  cudaGlobalDevice = -1;

  if (use_cuda) {

    pid_t pid = getpid();

    // Get device count; exit on failure
    //
    cuda_safe_call_mpi(cudaGetDeviceCount(&deviceCount), __FILE__, __LINE__,
		       myid, "cudaGetDevicecCount failure");

    // Query and assign my CUDA device
    //
    if (deviceCount>0) {

      int totalCount = std::max<int>(deviceCount, ngpus);

      // Get my local rank in sibling processes
      //
      int myCount = 0, curCount = 0;
      for (auto v : siblingList) {
	if (myid==v) myCount = curCount;
	curCount++;
      }
      
      // Allow GPU to be used by multiple MPI processes
      //
      if (myCount < totalCount) cudaGlobalDevice = myCount % deviceCount;
      
      // Set device; exit on failure
      //
      if (cudaGlobalDevice>=0) {

	cuda_safe_call_mpi(cudaSetDevice(cudaGlobalDevice), __FILE__, __LINE__,
			   myid, "cudaSetDevice failure");

	std::cout << "---- Process <" << pid << ">: "
		  << "setting CUDA device on Rank [" << myid
		  << "] on [" << processor_name << "] to ["
		  << cudaGlobalDevice << "/" << deviceCount << "]"
		  << std::endl;

      } else {
	
	std::cout << "---- Component <" << pid << ">: "
		  << "could not set CUDA device on Rank [" << myid
		  << "] on [" << processor_name << "] . . . "
		  << "this will cause a failure" << std::endl;
	
      }

    } else {
      std::ostringstream sout;
      sout << "[#" << myid << "] CUDA detected but deviceCount<=0!";
      throw GenericError(sout.str(), __FILE__, __LINE__, 1000, true);
    }

  }
#endif
}
