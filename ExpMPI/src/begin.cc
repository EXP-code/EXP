/*
  Get the whole thing going by:
  -----------------------------
  1) reading in initial phase space
  2) computing potential initially
  3) dumping first phase-space and log file
  4) write parameter file
*/

#include "expand.h"

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void begin_run(void)
{

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
  // Initialize phase-space components 
  //===================================

  comp.initialize();

  //===================================
  // Initialize external forces, if any
  //===================================

  external.initialize();
  
  //===============================
  // Compute initial accereration  
  //===============================

#ifdef DEBUG
  cout << "Process " << myid << ": about to compute potential [begin]\n";
#endif
  comp.compute_potential();
#ifdef DEBUG
  cout << "Process " << myid << ": potential computed [begin]\n";
#endif

  //==============================
  // Increment velocity of bodies 
  // to initialize leap-frog      
  //==============================

  init_velocity();

  //===================================
  // Initialize output routines
  //===================================

  output.initialize();
  output.Run(0);

  //======================
  // Write parameter file 
  //======================

  write_parm();

}
