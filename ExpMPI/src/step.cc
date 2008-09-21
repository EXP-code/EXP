/* 
  Call necessary routines to advance phase-space one step
*/

#include "expand.h"

#include <Timer.h>
static Timer timer_coef(true), timer_drift(true), timer_vel(true);
static Timer timer_pot(true), timer_adj(true);
static unsigned tskip = 1;
static bool timing = false;

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

void do_step(int n)
{
  // Turn on step timers or VERBOSE level 4 or greater
  //
  if (VERBOSE>3) timing = true;

  //========================
  // Advance using leapfrog 
  // algorithm:
  //
  // K_{1/2} D_1 K_{1/2}
  //========================

  comp.multistep_reset();

  if (comp.bad_values()) 
    cout << "Process " << myid 
	 << ": found BAD values before multistep" << endl;

  if (multistep) {
    
				// Sub step markers
    for (int M=0; M<=multistep; M++) stepL[M] = stepN[M] = 0;

    double dt = dtime/Mstep;	// Smallest time step

				// COM update
    if (timing) timer_vel.start();
    incr_com_velocity(0.5*dt);
    if (timing) timer_vel.stop();

				// March through all the substeps
				// of the hierarchy
    for (mstep=1; mstep<=Mstep; mstep++) {

      for (int M=0; M<=multistep; M++) {

	if (timing) timer_coef.start();
	if (stepN[M]<mstep) {
				// Time step at this level
				// 
	  double DT = dt*mintvl[M];

				// Advance velocity by 1/2 step:
				// First K_{1/2}
	  if (timing) timer_vel.start();
	  incr_velocity(0.5*DT);
	  if (timing) timer_vel.stop();

	  if (comp.bad_values()) 
	    cout << "Process " << myid 
		 << ": found BAD values after incr_vel M=" << M << endl;

				// Advance position by whole step:
				// D_1
	  if (timing) timer_drift.start();
	  incr_position(DT, M);
	  if (timing) timer_drift.stop();

	  if (comp.bad_values()) 
	    cout << "Process " << myid 
		 << ": found BAD values after incr_posl M=" << M << endl;

				// Do the swap
				//
	  comp.multistep_swap(M);

	  if (comp.bad_values()) 
	    cout << "Process " << myid 
		 << ": found BAD values after swap M=" << M << endl;

				// Update the markers
				// 
	  stepL[M]  =  stepN[M];
	  stepN[M] += mintvl[M];

				// Compute the coefficients
				// for this level
	  comp.compute_expansion(M);

	  if (comp.bad_values()) 
	    cout << "Process " << myid 
		 << ": found BAD values after expansion M=" << M << endl;

				// Reset the particle positions;
				// they will be drifted sycrhonously below
	  if (timing) timer_drift.start();
	  if (posnsync) incr_position(-DT, M);
	  if (timing) timer_drift.stop();

	  if (comp.bad_values()) 
	    cout << "Process " << myid 
		 << ": found BAD values after second incr pos M=" << M << endl;
	  /*
	  if (myid==0) {
	    if (mstep==1 && M==0) cout << endl;
	    cout << "Step=" << setw(3) << n
		 << " | Substep=" << setw(3) << mstep
		 << " | Mlevel=" << setw(3) << M
		 << " | Mfirst=" << setw(3) << mfirst[mstep-1]
		 << " | deltaM=" << setw(3) << mintvl[M]
		 << " | stepL=" << setw(3) << stepL[M]
		 << " | stepN=" << setw(3) << stepN[M]
		 << " | T=" << setw(18) << tnow
		 << " | DT=" << setw(18) << DT
		 << " | advT=" << setw(18) << tnow+DT
		 << endl;
	  }
	  */
	}
	if (timing) timer_coef.stop();

	if (comp.bad_values())  
	  cout << "Process " << myid
	       << ": found BAD values after multistep coef" << endl;
      }

      tnow += dt;		// Time at the end of the current step

				// Drift all positions by the substep size
      if (timing) timer_drift.start();
      if (posnsync) incr_position(dt, -1);
      //                              ^
      //                              |
      //           ignoring levels----/

				// COM update
      incr_com_position(dt);
      if (timing) timer_drift.stop();

				// Compute potential at active levels
      if (timing) timer_pot.start();
      comp.compute_potential(mfirst[mstep-1]);
      if (timing) timer_pot.stop();

      if (comp.bad_values())
	cout << "Process " << myid
	     << ": found BAD values after compute_potential" << endl;

				// For all active levels . . .
				// Advance velocity by 1/2 step:
				// Second K_{1/2}
      if (timing) timer_vel.start();
      for (int M=mfirst[mstep-1]; M<=multistep; M++) {
	incr_velocity(0.5*dt*mintvl[M], M);
	// For low-level debugging . . .
	/*
	if (myid==0) {
	  cout << "Step=" << setw(3) << n
	       << " | Substep=" << setw(3) << mstep
	       << " | Mlevel=" << setw(3) << M
	       << " | Mfirst=" << setw(3) << mfirst[mstep-1]
	       << " | velocity advanced by DT=" << 0.5*dt*mintvl[M]
	       << endl;
	}
	*/
      }
				// COM update
      incr_com_velocity(0.5*dt);
      if (timing) timer_vel.stop();

      if (timing) timer_adj.start();
      adjust_multistep_level(false);
      if (timing) timer_adj.stop();

      if (comp.bad_values())
	cout << "Process " << myid
	     << ": found BAD values after multistep advance" << endl;
    }

  } else {
				// Time at the end of the step
    tnow += dtime;
				// Velocity by 1/2 step
    if (timing) timer_vel.start();
    incr_velocity(0.5*dtime);
    incr_com_velocity(0.5*dtime);
    if (timing) timer_vel.stop();
				// Position by whole step
    if (timing) timer_drift.start();
    incr_position(dtime);
    incr_com_position(dtime);
    if (timing) timer_drift.start();
				// Compute acceleration
    if (timing) timer_pot.start();
    comp.compute_potential();
    if (timing) timer_pot.stop();
				// Velocity by 1/2 step
    if (timing) timer_vel.start();
    incr_velocity(0.5*dtime);
    incr_com_velocity(0.5*dtime);
    if (timing) timer_vel.stop();
  }

				// Write output
  output.Run(n);

				// Load balance
  comp.load_balance();

				// Timer output
  if (timing && this_step!=0 && (this_step % tskip) == 0) {
    if (myid==0) {
      cout << endl
	   << setw(70) << setfill('-') << '-' << endl
	   << setw(70) << left << "--- Timer info" << endl
	   << setw(70) << setfill('-') << '-' << endl << setfill(' ') << right;
      if (multistep) {
	cout << setw(20) << "Coefs: "
	     << setw(18) << 1.0e-6*timer_coef.getTime().getRealTime() << endl
	     << setw(20) << "Drift: "
	     << setw(18) << 1.0e-6*timer_drift.getTime().getRealTime() << endl
	     << setw(20) << "Velocity: "
	     << setw(18) << 1.0e-6*timer_vel.getTime().getRealTime() << endl
	     << setw(20) << "Force: "
	     << setw(18) << 1.0e-6*timer_pot.getTime().getRealTime() << endl
	     << setw(20) << "Adjust: "
	     << setw(18) << 1.0e-6*timer_adj.getTime().getRealTime() << endl;

      } else {
	cout << setw(20) << "Drift: "
	     << setw(18) << 1.0e-6*timer_drift.getTime().getRealTime() << endl
	     << setw(20) << "Velocity: "
	     << setw(18) << 1.0e-6*timer_vel.getTime().getRealTime() << endl
	     << setw(20) << "Force: "
	     << setw(18) << 1.0e-6*timer_pot.getTime().getRealTime() << endl;
      }
      cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
    }

    // DEBUG
    //
    if (VERBOSE>4) {
      
      vector<int> levpop(multistep+1, 0), levtot(multistep+1, 0);
      for (list<Component*>::iterator cc=comp.components.begin(); 
	   cc != comp.components.end(); cc++) {
	for (int n=0; n<=multistep; n++) levpop[n] += (*cc)->levlist[n].size();
      }

      MPI_Reduce(&levpop[0], &levtot[0], multistep+1, MPI_INT, MPI_SUM, 0,
		 MPI_COMM_WORLD);

      if (myid==0) {
	cout << endl
	     << setw(70) << setfill('-') << '-' << endl
	     << setw(70) << left << "--- Level info" << endl
	     << setw(70) << setfill('-') << '-' << endl 
	     << setfill(' ') << right;
      }
      for (int n=0; n<numprocs; n++) {
	if (myid==n) {
	  cout << setw(4) << myid << ": ";
	  for (int m=0; m<=multistep; m++) cout << setw(8) << levpop[m];
	  cout << endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
      }

      if (myid==0) {
	cout << endl;
	cout << setw(4) << "T" << ": ";
	for (int m=0; m<=multistep; m++) cout << setw(8) << levtot[m];
	cout << endl;
	cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
      }
    }
    //
    // END DEBUG

    timer_coef.reset();
    timer_drift.reset();
    timer_vel.reset();
    timer_pot.reset();
    timer_adj.reset();
  }
}
