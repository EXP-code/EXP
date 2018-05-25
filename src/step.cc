/**
  Call necessary routines to advance phase-space one step
*/

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

#if HAVE_LIBCUDA==1
#include <NVTX.H>
#endif

static Timer timer_coef(true), timer_drift(true), timer_vel(true);
static Timer timer_pot (true), timer_adj  (true);

static unsigned tskip = 1;
static bool timing = false;

inline void check_bad(const char *msg)
{
#ifdef CHK_BADV
  if (comp.bad_values()) 
    cout << "Process " << myid 
	 << ": found BAD values " << msg << endl;
#endif
}

inline void check_bad(const char *msg, int v)
{
#ifdef CHK_BADV
  if (comp.bad_values()) 
    cout << "Process " << myid 
	 << ": found BAD values " << msg << ", M=" << v << endl;
#endif
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

  comp->multistep_reset();

  check_bad("before multistep");

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
				// particles that move on this step
				// (the "active" particles)
				//
      for (int M=mfirst[mstep]; M<=multistep; M++) {
	nvTracerPtr tPtr2;

				// The timestep at level M
	double DT = dt*mintvl[M];
	
				// Advance velocity by 1/2 step:
				// First K_{1/2}
	if (cuda_prof) tPtr2 = nvTracerPtr(new nvTracer("Velocity kick [1]"));
	if (timing) timer_vel.start();
	incr_velocity(0.5*DT, M);
#ifdef CHK_STEP
	vel_check[M] += 0.5*DT;
#endif
	if (timing) timer_vel.stop();

	check_bad("after incr_vel", M);

				// Advance position by the whole time
				// step at this level: D_1
				//
	if (cuda_prof) tPtr2 = nvTracerPtr(new nvTracer("Drift"));
	if (timing) timer_drift.start();
	incr_position(DT, M);
#ifdef CHK_STEP
	pos_check[M] += DT;
#endif
	if (timing) timer_drift.stop();

	check_bad("after incr_pos", M);

				// Now, compute the coefficients for
				// this level
				//
	if (cuda_prof) tPtr2 = nvTracerPtr(new nvTracer("Expansion"));
	if (timing) timer_coef.start();
	comp->compute_expansion(M);
	if (timing) timer_coef.stop();
      }

      nvTracerPtr tPtr1;

      
      double tlast = tnow;	// Time before current step
      tnow += dt;		// Time at the end of the current step

				// COM update:
				// Position drift
      if (timing) timer_drift.start();
      incr_com_position(dt);
      if (timing) timer_drift.stop();

				// Compute potential for all the
				// particles active at this step
      if (cuda_prof) tPtr1 = nvTracerPtr(new nvTracer("Potential"));
      if (timing) timer_pot.start();
      comp->compute_potential(mfirst[mstep]);
      if (timing) timer_pot.stop();

      check_bad("after compute_potential");

				// For all active levels . . .
				// Advance velocity by 1/2 step to
				// bring the velocity in sync: 
				// Second K_{1/2}
      if (cuda_prof) tPtr1 = nvTracerPtr(new nvTracer("Velocity kick [2]"));
      if (timing) timer_vel.start();
      for (int M=mfirst[mstep]; M<=multistep; M++) {
	incr_velocity(0.5*dt*mintvl[M], M);
#ifdef CHK_STEP
	vel_check[M] += 0.5*dt*mintvl[M];
#endif
      }
      if (timing) timer_vel.stop();

      if (cuda_prof) tPtr1 = nvTracerPtr(new nvTracer("Adjust multistep"));
      if (timing) timer_adj.start();
      adjust_multistep_level(false);
      if (mstep==0) { // Print the level lists
	comp->print_level_lists(tlast);
      }
      if (timing) timer_adj.stop();

      check_bad("after multistep advance");

				// DEBUG
#ifdef DEBUG
      comp->multistep_debug();
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
    nvTracerPtr tPtr1;

				// Time at the end of the step
    tnow += dtime;
				// Velocity by 1/2 step
    if (cuda_prof) tPtr1 = nvTracerPtr(new nvTracer("Velocity kick [1]"));
    if (timing) timer_vel.start();
    incr_velocity(0.5*dtime);
    incr_com_velocity(0.5*dtime);
    if (timing) timer_vel.stop();
				// Position by whole step
    if (cuda_prof) tPtr1 = nvTracerPtr(new nvTracer("Drift"));
    if (timing) timer_drift.start();
    incr_position(dtime);
    incr_com_position(dtime);
    if (timing) timer_drift.start();
				// Compute acceleration
    if (cuda_prof) tPtr1 = nvTracerPtr(new nvTracer("Potential"));
    if (timing) timer_pot.start();
    comp->compute_potential();
    if (timing) timer_pot.stop();
				// Velocity by 1/2 step
    if (cuda_prof) tPtr1 = nvTracerPtr(new nvTracer("Velocity kick [2]"));
    if (timing) timer_vel.start();
    incr_velocity(0.5*dtime);
    incr_com_velocity(0.5*dtime);
    if (timing) timer_vel.stop();
  }

  nvTracerPtr tPtr;
				// Write output
  if (cuda_prof) tPtr = nvTracerPtr(new nvTracer("Data output"));
  output->Run(n);

				// Summarize processor particle load
  comp->report_numbers();

				// Load balance
  if (cuda_prof) tPtr = nvTracerPtr(new nvTracer("Load balance"));
  comp->load_balance();

				// Timer output
  if (timing && this_step!=0 && (this_step % tskip) == 0) {
    if (myid==0) {
      cout << endl
	   << setw(70) << setfill('-') << '-' << endl
	   << setw(70) << left << "--- Timer info" << endl
	   << setw(70) << setfill('-') << '-' << endl << setfill(' ') << right;
      if (multistep) {
	cout << setw(20) << "Coefs: "
	     << setw(18) << timer_coef.getTime()() << endl
	     << setw(20) << "Drift: "
	     << setw(18) << timer_drift.getTime()() << endl
	     << setw(20) << "Velocity: "
	     << setw(18) << timer_vel.getTime()() << endl
	     << setw(20) << "Force: "
	     << setw(18) << timer_pot.getTime()() << endl
	     << setw(20) << "Adjust: "
	     << setw(18) << timer_adj.getTime()() << endl;
      }
      cout << setw(70) << setfill('-') << '-' << endl << setfill(' ');
    }

    //
    // DEBUG
    //
    if (VERBOSE>4) {
      
      vector<int> levpop(multistep+1, 0), levtot(multistep+1, 0);
      for (auto c : comp->components) {
	for (int n=0; n<=multistep; n++) levpop[n] += c->levlist[n].size();
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
