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

#include <NVTX.H>

static Timer timer_coef, timer_drift, timer_vel;
static Timer timer_pot , timer_adj  , timer_tot;
static Timer timer_cuda;

static unsigned tskip = 1;

inline void check_bad(const char *msg)
{
#ifdef CHK_BADV
  if (comp.bad_values()) 
    std::cout << "Process " << myid 
	      << ": found BAD values " << msg << std::endl;
#endif
}

inline void check_bad(const char *msg, int v)
{
#ifdef CHK_BADV
  if (comp.bad_values()) 
    std::cout << "Process " << myid 
	      << ": found BAD values " << msg << ", M=" << v << std::endl;
#endif
}

void do_step(int n)
{
#ifdef USE_GPTL
  GPTLstart("dostep");
#endif

  // Turn on step timers or VERBOSE level 4 or greater
  //
  if (VERBOSE>3) step_timing = true;

  //========================
  // Advance using leapfrog 
  // algorithm:
  //
  // K_{1/2} D_1 K_{1/2}
  //========================

  comp->multistep_reset();

  check_bad("before multistep");

  if (step_timing) timer_tot.start();

  // set up CUDA tracer
  nvTracerPtr tPtr;

  // BEG: multistep>0 block
  if (multistep) {
    
    double dt = dtime/Mstep;	// Smallest time step

				// COM update:
				// First velocity half-kick
    if (step_timing) timer_vel.start();
    incr_com_velocity(0.5*dtime); 
    if (step_timing) timer_vel.stop();

#ifdef CHK_STEP
    vector<double> pos_check(multistep+1);
    vector<double> vel_check(multistep+1);
#endif
    
				// March through all the substeps
				// of the hierarchy
    for (mstep=0; mstep<Mstep; mstep++) {

				// Write multistep output
      if (cuda_prof) tPtr = nvTracerPtr(new nvTracer("Data output"));
      output->Run(n, mstep);

				// Compute next coefficients for
				// particles that move on this step
				// (the "active" particles)
				//
      for (int M=mfirst[mstep]; M<=multistep; M++) {
				// The timestep at level M
	double DT = dt*mintvl[M];
	
				// Advance velocity by 1/2 step:
				// First K_{1/2}
	nvTracerPtr tPtr2;
	if (cuda_prof) {
	  tPtr2 = nvTracerPtr(new nvTracer("Velocity kick [1]"));
	}
	if (step_timing) timer_vel.start();
	incr_velocity(0.5*DT, M);
#ifdef CHK_STEP
	vel_check[M] += 0.5*DT;
#endif
	if (step_timing) timer_vel.stop();

	check_bad("after incr_vel", M);

				// Advance position by the whole time
				// step at this level: D_1
				//
	if (cuda_prof) {
	  tPtr2.reset();
	  tPtr2 = nvTracerPtr(new nvTracer("Drift"));
	}
	if (step_timing) timer_drift.start();
	incr_position(DT, M);
#ifdef CHK_STEP
	pos_check[M] += DT;
#endif
	if (step_timing) timer_drift.stop();

	check_bad("after incr_pos", M);

				// Now, compute the coefficients for
				// this level
				//

	if (cuda_prof) {
	  tPtr2.reset();
	  tPtr2 = nvTracerPtr(new nvTracer("Expansion"));
	}
	if (step_timing) timer_coef.start();
	comp->compute_expansion(M);
	if (step_timing) timer_coef.stop();
      }
      
				// COM update:
				// Position drift
      if (step_timing) timer_drift.start();
      incr_com_position(dt);
      if (step_timing) timer_drift.stop();

				// Compute potential for all the
				// particles active at this step
      nvTracerPtr tPtr1;
      if (cuda_prof) tPtr1 = nvTracerPtr(new nvTracer("Potential"));
      if (step_timing) timer_pot.start();
      {
	double tlast = tnow;	// Time before current step
				// Time at the end of the drift
	tstp = tnow + dt*mintvl[mfirst[mstep]];
				// Compute potential at drifted time
	comp->compute_potential(mfirst[mstep]);
	tstp  = tlast;		// Restore time to beginning of step
      }
      if (step_timing) timer_pot.stop();

      check_bad("after compute_potential");

				// For all active levels . . .
				// Advance velocity by 1/2 step to
				// bring the velocity in sync: 
				// Second K_{1/2}

      if (cuda_prof) {
	tPtr1.reset();
	tPtr1 = nvTracerPtr(new nvTracer("Velocity kick [2]"));
      }
      if (step_timing) timer_vel.start();
      for (int M=mfirst[mstep]; M<=multistep; M++) {
	incr_velocity(0.5*dt*mintvl[M], M);
#ifdef CHK_STEP
	vel_check[M] += 0.5*dt*mintvl[M];
#endif
      }
      if (step_timing) timer_vel.stop();

      check_bad("after multistep advance");
				// DEBUG
#ifdef DEBUG
      comp->multistep_debug();
#endif
      if (mstep==0) {
				// Do particles at top level
	if (step_timing) timer_adj.start();
	adjust_multistep_level(true);
	if (step_timing) timer_adj.stop();
	
				// Print the level lists
	comp->print_level_lists(tnow);
      } else {
				// Do particles at lower levels
	if (step_timing) timer_adj.start();
	adjust_multistep_level(false);
	if (step_timing) timer_adj.stop();

      }

      tnow += dt;		// Next substep
      tstp  = tnow;		// Sync force time
    }
    // END: mstep loop

    // Write output
    if (cuda_prof) tPtr = nvTracerPtr(new nvTracer("Data output"));
    output->Run(n);

    if (cuda_prof) {
      tPtr = nvTracerPtr(new nvTracer("Adjust multistep"));
    }

				// COM update:
				// Second velocity half-kick
    if (step_timing) timer_vel.start();
    incr_com_velocity(0.5*dtime);
    if (step_timing) timer_vel.stop();

#ifdef CHK_STEP
				// Check steps
    if (myid==0) {		// 
      bool chk = true;
      for (int M=0; M<=multistep; M++) {
	if (fabs(dtime - pos_check[M]) > 1.0e-8*dtime) {
	  cerr << "Pos step error[" << M << "]: T=" << tnow << " found="
	       << pos_check[M] << ", expected=" << dtime << std::endl;
	  chk = false;
	}
	if (fabs(dtime - vel_check[M]) > 1.0e-8*dtime) {
	  cerr << "Vel step error[" << M << "]: T=" << tnow << " found="
	       << pos_check[M] << ", expected=" << dtime << std::endl;
	  chk = false;
	}
      }
      if (chk) cerr << "Incremental steps OK at T=" << tnow << std::endl;
    }
#endif

  }
  // END: multistep>0 block
  // BEG: multistep=0 block
  else {
				// Time at the end of the step
    tnow += dtime;
    tstp  = tnow;
				// Velocity by 1/2 step
    nvTracerPtr tPtr1;
    if (cuda_prof) tPtr1 = nvTracerPtr(new nvTracer("Velocity kick [1]"));
    if (step_timing) timer_vel.start();
    incr_velocity(0.5*dtime);
    incr_com_velocity(0.5*dtime);
    if (step_timing) timer_vel.stop();
				// Position by whole step
    if (cuda_prof) {
      tPtr1.reset();
      tPtr1 = nvTracerPtr(new nvTracer("Drift"));
    }
    if (step_timing) timer_drift.start();
    incr_position(dtime);
    incr_com_position(dtime);
    if (step_timing) timer_drift.stop();

				// Compute coefficients
    if (step_timing) timer_coef.start();
    comp->compute_expansion(0);
    if (step_timing) timer_coef.stop();

				// Compute acceleration
    if (cuda_prof) {
      tPtr1.reset();
      tPtr1 = nvTracerPtr(new nvTracer("Potential"));
    }
    if (step_timing) timer_pot.start();
    comp->compute_potential();
    if (step_timing) timer_pot.stop();
				// Velocity by 1/2 step
    if (cuda_prof) {
      tPtr1.reset();
      tPtr1 = nvTracerPtr(new nvTracer("Velocity kick [2]"));
    }
    if (step_timing) timer_vel.start();
    incr_velocity(0.5*dtime);
    incr_com_velocity(0.5*dtime);
    if (step_timing) timer_vel.stop();

                                 // Write output
    nvTracerPtr tPtr;
    if (cuda_prof) tPtr = nvTracerPtr(new nvTracer("Data output"));
    output->Run(n);

  }
  // END: multistep=0 block

  if (step_timing) timer_tot.stop();

				// Write output
  //nvTracerPtr tPtr;
  //if (cuda_prof) tPtr = nvTracerPtr(new nvTracer("Data output"));
  //output->Run(n);

				// Summarize processor particle load
  comp->report_numbers();

				// Load balance
  if (cuda_prof) {
    tPtr.reset();
    tPtr = nvTracerPtr(new nvTracer("Load balance"));
  }
  comp->load_balance();

				// Timer output
  if (step_timing && this_step!=0 && (this_step % tskip) == 0) {
    if (myid==0) {
      auto totalT = timer_tot.getTime();
      std::cout << std::endl
		<< std::setw(70) << std::setfill('-') << '-' << std::endl
		<< std::setw(70) << left << "--- Timer info" << std::endl
		<< std::setw(70) << std::setfill('-') << '-' << std::endl
		<< std::setfill(' ') << std::right
		<< std::setw(20) << "Drift: "
		<< std::setw(18) << timer_drift.getTime()
		<< std::setw(18) << timer_drift.getTime()/totalT << std::endl
		<< std::setw(20) << "Velocity: "
		<< std::setw(18) << timer_vel.getTime()
		<< std::setw(18) << timer_vel.getTime()/totalT << std::endl
		<< std::setw(20) << "Force: "
		<< std::setw(18) << timer_pot.getTime()
		<< std::setw(18) << timer_pot.getTime()/totalT << std::endl
		<< std::setw(20) << "Coefs: "
		<< std::setw(18) << timer_coef.getTime()
		<< std::setw(18) << timer_coef.getTime()/totalT << std::endl;
      if (multistep)
	std::cout << std::setw(20) << "Adjust: "
		  << std::setw(18) << timer_adj.getTime()
		  << std::setw(18) << timer_adj.getTime()/totalT << std::endl;
      if (use_cuda)
	std::cout << std::setw(20) << "Cuda copy: "
		  << std::setw(18) << comp->timer_cuda.getTime()
		  << std::setw(18) << comp->timer_cuda.getTime()/totalT << std::endl;
      std::cout << std::setw(20) << "Total: "
		<< std::setw(18) << timer_tot.getTime()
		<< std::setw(18) << 1.0 << std::endl
		<< std::setw(70) << std::setfill('-') << '-' << std::endl
		<< std::setfill(' ');
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
	std::cout << std::endl
		  << std::setw(70) << std::setfill('-') << '-' << std::endl
		  << std::setw(70) << left << "--- Level info" << std::endl
		  << std::setw(70) << std::setfill('-') << '-' << std::endl 
		  << std::setfill(' ') << right;
      }
      for (int n=0; n<numprocs; n++) {
	if (myid==n) {
	  std::cout << std::setw(4) << myid << ": ";
	  for (int m=0; m<=multistep; m++) std::cout << std::setw(8) << levpop[m];
	  std::cout << std::endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);
      }

      if (myid==0) {
	std::cout << std::endl;
	std::cout << std::setw(4) << "T" << ": ";
	for (int m=0; m<=multistep; m++) std::cout << std::setw(8) << levtot[m];
	std::cout << std::endl;
	std::cout << std::setw(70) << std::setfill('-') << '-' << std::endl
		  << std::setfill(' ');
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
    timer_tot  .reset();
    if (use_cuda) comp->timer_cuda.reset();
  }

#ifdef USE_GPTL
  GPTLstop("dostep");
#endif
}
