/*
  Get the whole thing going by:
  -----------------------------
  1) reading in initial phase space
  2) computing potential initially
  3) dumping first phase-space and log file
  4) write parameter file
*/

#include <expand.h>
#include <ExternalCollection.H>
#include <OutputContainer.H>

void begin_run(void)
{

  //===================================
  // Make the instance containers
  //===================================
  
  comp     = new ComponentContainer;
  external = new ExternalCollection;
  output   = new OutputContainer;

  //===================================
  // Make the Vector class mutex
  //===================================
  
  int errcode = pthread_mutex_init(&mem_lock, NULL);
  if (errcode) {
    cerr << "Process " << myid 
	 << ": failed to make the Vector class memory lock, errcode=" 
	 << errcode << endl;
    MPI_Abort(MPI_COMM_WORLD, 115);
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
  
  if (multistep) {
    comp->multistep_reset();

    comp->compute_expansion(0);
    //                      ^
    //                      |
    //                      |
    comp->compute_potential(0);
    //                      ^
    //                      |
    //   All time levels----/

    //==============================
    // Initialize multistep levels
    //==============================
    adjust_multistep_level(true);
    //                     ^
    //                     |
    // Do all particles----+
  }

  // Compute coefficients . . . 
  //
  if (multistep) comp->multistep_reset();

  comp->compute_expansion(0);
  comp->compute_potential(0);

  initializing = false;

  if (multistep) comp->multistep_reset();

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
