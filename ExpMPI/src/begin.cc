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

  comp.compute_potential();

  //===================================
  // Initialize output routines
  //===================================

  output.initialize();
  output.Run(0);

  //==============================
  // Increment velocity of bodies 
  // to initialize leap-frog      
  //==============================

  init_velocity();

  //======================
  // Write parameter file 
  //======================

  write_parm();

}
