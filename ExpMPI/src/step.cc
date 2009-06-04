/* 
  Call necessary routines to advance phase-space one step
*/

#ifdef RCSID
static char rcsid[] = "$Id$";
#endif

// Uncomment for time step debugging
//
// #define CHK_STEP

// Uncomment for bad value checking
//
// #define CHK_BADV

#include <expand.h>
#include <OutputContainer.H>

// Substep timing
//
#include <Timer.h>

#ifdef USE_GPTL
#include <gptl.h>
#endif

static Timer timer_coef(true), timer_drift(true), timer_vel(true);
static Timer timer_pot(true), timer_adj(true);
static unsigned tskip = 1;
static bool timing = false;

void check_bad(const char *msg)
{
  if (comp.bad_values()) 
    cout << "Process " << myid 
	 << ": found BAD values " << msg << endl;
}

void check_bad(const char *msg, int v)
{
  if (comp.bad_values()) 
    cout << "Process " << myid 
	 << ": found BAD values " << msg << ", M=" << v << endl;
}

void do_step(int n)
{
#ifdef USE_GPTL
  GPTLstart("dostep");
#endif


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

#ifdef CHK_BADV
  check_bad("before multistep");
#endif

  if (multistep) {
    
    double dt = dtime/Mstep;	// Smallest time step

				// COM update:
				// First velocity half-kick
    if (timing) timer_vel.start();
    incr_com_velocity(0.5*dtime);
    if (timing) timer_vel.stop();

#ifdef CHK_STEP
    vector<double> pos_check(multistep+1);
    vector<double> vel_check(multistep+1);
#endif
    
				// March through all the substeps
				// of the hierarchy
    for (mstep=0; mstep<Mstep; mstep++) {

				// Compute next coefficients for
				// active steps
      for (int M=mfirst[mstep]; M<=multistep; M++) {

				// 
	double DT = dt*mintvl[M];
	
				// Advance velocity by 1/2 step:
				// First K_{1/2}
	if (timing) timer_vel.start();
	incr_velocity(0.5*DT, M);
#ifdef CHK_STEP
	vel_check[M] += 0.5*DT;
#endif
	if (timing) timer_vel.stop();

#ifdef CHK_BADV
	check_bad("after incr_vel", M);
#endif

				// Advance position by whole step:
				// D_1
	if (timing) timer_drift.start();
	incr_position(DT, M);
#ifdef CHK_STEP
	pos_check[M] += DT;
#endif
	if (timing) timer_drift.stop();

#ifdef CHK_BADV
	check_bad("after incr_pos", M);
#endif
				// Compute the coefficients
				// for this level
	if (timing) timer_coef.start();
	comp.compute_expansion(M);
	if (timing) timer_coef.stop();

#ifdef CHK_BADV
	check_bad("after expansion", M);
#endif	
				// Reset the particle positions by a whole step
				// they will be drifted sycrhonously below
	if (timing) timer_drift.start();
	if (posnsync) {
	  incr_position(-DT, M);
#ifdef CHK_STEP
	  pos_check[M] -= DT;
#endif
	}
	if (timing) timer_drift.stop();

#ifdef CHK_BADV
	check_bad("after second incr pos", M);
#endif
      }

      
      double tlast = tnow;	// Time before current step
      tnow += dt;		// Time at the end of the current step

				// Drift all positions by the substep size
      if (timing) timer_drift.start();
#ifdef CHK_STEP
      if (posnsync) {
	for (int M=0; M<=multistep; M++) pos_check[M] += dt;
      }
#endif
      if (posnsync) incr_position(dt, -1);
      //                              ^
      //                              |
      //           ignoring levels----/

				// COM update:
				// Position drift
      incr_com_position(dt);
      if (timing) timer_drift.stop();

				// Compute potential at active levels
				// 
      if (timing) timer_pot.start();
      comp.compute_potential(mfirst[mstep]);
      if (timing) timer_pot.stop();

#ifdef CHK_BADV
      check_bad("after compute_potential");
#endif

				// For all active levels . . .
				// Advance velocity by 1/2 step:
				// Second K_{1/2}
      if (timing) timer_vel.start();
      for (int M=mfirst[mstep]; M<=multistep; M++) {
	incr_velocity(0.5*dt*mintvl[M], M);
#ifdef CHK_STEP
	vel_check[M] += 0.5*dt*mintvl[M];
#endif
      }
      if (timing) timer_vel.stop();

      if (timing) timer_adj.start();
      adjust_multistep_level(false);
      if (mstep==0) {		// Print the level lists
	comp.print_level_lists(tlast);
      }
      if (timing) timer_adj.stop();

#ifdef CHK_BADV
      check_bad("after multistep advance");
#endif

				// DEBUG
#ifdef DEBUG
      comp.multistep_debug();
#endif
    }

				// COM update:
				// Second velocity half-kick
    if (timing) timer_vel.start();
    incr_com_velocity(0.5*dtime);
    if (timing) timer_vel.stop();

#ifdef CHK_STEP
				// Check steps
    if (myid==0) {		// 
      bool chk = true;
      for (int M=0; M<=multistep; M++) {
	if (fabs(dtime - pos_check[M]) > 1.0e-8*dtime) {
	  cerr << "Pos step error[" << M << "]: T=" << tnow << " found="
	       << pos_check[M] << ", expected=" << dtime << endl;
	  chk = false;
	}
	if (fabs(dtime - vel_check[M]) > 1.0e-8*dtime) {
	  cerr << "Vel step error[" << M << "]: T=" << tnow << " found="
	       << pos_check[M] << ", expected=" << dtime << endl;
	  chk = false;
	}
      }
      if (chk) cerr << "Incremental steps OK at T=" << tnow << endl;
    }
#endif

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
      }
      cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
    }

    //
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
    //

    timer_coef .reset();
    timer_drift.reset();
    timer_vel  .reset();
    timer_pot  .reset();
    timer_adj  .reset();
  }

#ifdef USE_GPTL
  GPTLstop("dostep");
#endif
}
