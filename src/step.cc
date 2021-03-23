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
  //
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
    
    // March through all the substeps of the hierarchy
    //
    for (mstep=0; mstep<Mstep; mstep++) {

				// Write multistep output
      if (cuda_prof) tPtr = nvTracerPtr(new nvTracer("Data output"));
      output->Run(n, mstep);

      // Compute next coefficients for particles that move on this
      // step (the "active" particles)
      //
      for (int M=mfirst[mstep]; M<=multistep; M++) {
				// The timestep at level M
	double DT = dt*mintvl[M];
	
	// Advance velocity by 1/2 step for active particles: First
	// K_{1/2}
	//
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

	// Advance position by the whole time step at this level: D_1
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

	// Now, compute the coefficients for this level at the
	// advanced position in preparation for the next kick
	//
	if (cuda_prof) {
	  tPtr2.reset();
	  tPtr2 = nvTracerPtr(new nvTracer("Expansion"));
	}
	if (step_timing) timer_coef.start();
	comp->compute_expansion(M);
	if (step_timing) timer_coef.stop();
      }
      
      tnow += dt;		// Time at the end of the current step

      // COM update: Position drift
      //
      if (step_timing) timer_drift.start();
      incr_com_position(dt);
      if (step_timing) timer_drift.stop();

      // Compute potential for all the particles active at this step
      //
      nvTracerPtr tPtr1;
      if (cuda_prof) tPtr1 = nvTracerPtr(new nvTracer("Potential"));
      if (step_timing) timer_pot.start();
      comp->compute_potential(mfirst[mstep]);
      if (step_timing) timer_pot.stop();

      check_bad("after compute_potential");

      // Second K_{1/2}
      // 
      // Advance velocity by 1/2 step for active particles at the next
      // sub step.  Inactive particles' velocities remain unsynced
      // otherwise because the coefficient values will not be
      // available all levels for their drifted positions until later
      // mstep.
      //
      if (cuda_prof) {
	tPtr1.reset();
	tPtr1 = nvTracerPtr(new nvTracer("Velocity kick [2]"));
      }

      if (step_timing) timer_vel.start();
      for (int M=mfirst[mstep+1]; M<=multistep; M++) {
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
    }
    // END: mstep loop

    // Write output
    if (cuda_prof) tPtr = nvTracerPtr(new nvTracer("Data output"));
    output->Run(n);

    if (cuda_prof) {
      tPtr = nvTracerPtr(new nvTracer("Adjust multistep"));
    }

    // COM update: second velocity half-kick
    //
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
      std::ostringstream sout;
      sout << "--- Timer info [T=" << tnow << "] ";
      std::cout << std::endl
		<< std::setw(70) << std::setfill('-') << '-' << std::endl
		<< std::setw(70) << left << sout.str()       << std::endl
		<< std::setw(70) << std::setfill('-') << '-' << std::endl
		<< std::setfill(' ') << std::right
		<< std::setw(20) << "Drift: " << std::scientific
		<< std::setw(18) << timer_drift.getTime() << std::fixed
		<< std::setw(18) << timer_drift.getTime()/totalT << std::endl
		<< std::setw(20) << "Velocity: " << std::scientific
		<< std::setw(18) << timer_vel.getTime() << std::fixed
		<< std::setw(18) << timer_vel.getTime()/totalT << std::endl
		<< std::setw(20) << "Force: " << std::scientific
		<< std::setw(18) << timer_pot.getTime() << std::fixed
		<< std::setw(18) << timer_pot.getTime()/totalT << std::endl
		<< std::setw(20) << "Coefs: " << std::scientific
		<< std::setw(18) << timer_coef.getTime() << std::fixed
		<< std::setw(18) << timer_coef.getTime()/totalT << std::endl;
      if (multistep)
	std::cout << std::setw(20) << "Adjust: " << std::scientific
		  << std::setw(18) << timer_adj.getTime() << std::fixed
		  << std::setw(18) << timer_adj.getTime()/totalT << std::endl;
      if (use_cuda) {
	std::cout << std::setw(20) << "Cuda copy: " << std::scientific
		  << std::setw(18) << comp->timer_cuda.getTime() << std::fixed
		  << std::setw(18) << comp->timer_cuda.getTime()/totalT << std::endl;
	std::cout << std::setw(20) << "Orient: " << std::scientific
		  << std::setw(18) << comp->timer_orient.getTime() << std::fixed
		  << std::setw(18) << comp->timer_orient.getTime()/totalT << std::endl;
      }
      std::cout << std::setw(20) << "Total: " << std::scientific
		<< std::setw(18) << timer_tot.getTime() << std::fixed
		<< std::setw(18) << 1.0 << std::endl << std::scientific
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
    if (use_cuda) comp->timer_orient.reset();
  }

#ifdef USE_GPTL
  GPTLstop("dostep");
#endif
}
