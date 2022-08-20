/**
  Call necessary routines to advance phase-space one step
*/

#include <memory>

// Uncomment for time step debugging
//
// #define CHK_STEP

// Uncomment for bad value checking
//
// #define CHK_BADV

#include <expand.H>
#include <OutputContainer.H>

// Substep timing
//
#include <Timer.H>

#ifdef USE_GPTL
#include <gptl.h>
#endif

#include <NVTX.H>

static Timer timer_coef, timer_drift, timer_vel, timer_out, timer_lev;
static Timer timer_pot , timer_adj  , timer_tot, timer_bal, timer_rpt;

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

// Hide all of the format manipulators here
std::string sForm(const std::string& label,
		  const double time, const double total)
{
  double pcnt = 0.0;
  if (total>0) pcnt = time/total*100.0;
  
  std::ostringstream str;
  str << std::setw(20) << std::left     << label << ": "
      << std::setw(18) << std::left     << time  << " ["
      << std::setw( 6) << std::right    << std::setprecision(2)
      << std::fixed    << pcnt << " %]" << std::endl;

  return str.str();
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

  // Start the total step timer
  //
  if (step_timing) timer_tot.start();

  // set up CUDA tracer
  //
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

      mdrft = mstep;		// Assign velocity position

				// Write multistep output
      if (step_timing) timer_out.start();
      if (cuda_prof) tPtr = std::make_shared<nvTracer>("Data output");
      output->Run(n, mstep);
      if (step_timing) timer_out.stop();

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
	  tPtr2 = std::make_shared<nvTracer>("Velocity kick [1]");
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
	  tPtr2 = std::make_shared<nvTracer>("Drift");
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
	  tPtr2 = std::make_shared<nvTracer>("Expansion");
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
      if (cuda_prof) tPtr1 = std::make_shared<nvTracer>("Potential");
      if (step_timing) timer_pot.start();
      mdrft = mstep + 1;	// Drifted position in multistep array
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
	tPtr1 = std::make_shared<nvTracer>("Velocity kick [2]");
      }

      if (step_timing) timer_vel.start();
      for (int M=mfirst[mdrft]; M<=multistep; M++) {
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
				// Adjust particle time-step levels
      if (step_timing) timer_adj.start();
      adjust_multistep_level(false);
      if (step_timing) timer_adj.stop();
      
      // Print the level lists
      if (step_timing) timer_lev.start();
      if (mdrft==Mstep) comp->print_level_lists(tnow);
      if (step_timing) timer_lev.stop();
    }
    // END: mstep loop

    // Write output
    if (step_timing) timer_out.start();
    if (cuda_prof) tPtr = std::make_shared<nvTracer>("Data output");
    output->Run(n);
    if (step_timing) timer_out.stop();

    if (cuda_prof) {
      tPtr = std::make_shared<nvTracer>("Adjust multistep");
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
    if (cuda_prof) tPtr1 = std::make_shared<nvTracer>("Velocity kick [1]");
    if (step_timing) timer_vel.start();
    incr_velocity(0.5*dtime);
    incr_com_velocity(0.5*dtime);
    if (step_timing) timer_vel.stop();
				// Position by whole step
    if (cuda_prof) {
      tPtr1.reset();
      tPtr1 = std::make_shared<nvTracer>("Drift");
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
      tPtr1 = std::make_shared<nvTracer>("Potential");
    }
    if (step_timing) timer_pot.start();
    comp->compute_potential();
    if (step_timing) timer_pot.stop();
				// Velocity by 1/2 step
    if (cuda_prof) {
      tPtr1.reset();
      tPtr1 = std::make_shared<nvTracer>("Velocity kick [2]");
    }
    if (step_timing) timer_vel.start();
    incr_velocity(0.5*dtime);
    incr_com_velocity(0.5*dtime);
    if (step_timing) timer_vel.stop();

                                 // Write output
    if (step_timing) timer_out.start();
    nvTracerPtr tPtr;
    if (cuda_prof) tPtr = std::make_shared<nvTracer>("Data output");
    output->Run(n);
    if (step_timing) timer_out.stop();

  }
  // END: multistep=0 block

				// Summarize processor particle load

  if (step_timing) timer_rpt.start();
  comp->report_numbers();
  if (step_timing) timer_rpt.stop();

				// Load balance
  if (cuda_prof) {
    tPtr.reset();
    tPtr = std::make_shared<nvTracer>("Load balance");
  }
  if (step_timing) timer_bal.start();
  comp->load_balance();
  if (step_timing) timer_bal.stop();

				// Stop the total step timer
  if (step_timing) timer_tot.stop();

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
		<< sForm("Drift",    timer_drift.getTime(), totalT)
		<< sForm("Velocity", timer_vel.getTime(),   totalT)
		<< sForm("Force",    timer_pot.getTime(),   totalT)
		<< sForm("Coefs",    timer_coef.getTime(),  totalT)
		<< sForm("Output",   timer_out.getTime(),   totalT)
		<< sForm("Levels",   timer_lev.getTime(),   totalT)
		<< sForm("Report",   timer_rpt.getTime(),   totalT)
		<< sForm("Balance",  timer_bal.getTime(),   totalT);
      if (multistep)
	std::cout << sForm("Adjust", timer_adj.getTime(),   totalT);
      if (use_cuda) {
	std::cout << sForm("Cuda copy", comp->timer_cuda.getTime(),   totalT)
		  << sForm("Orient",    comp->timer_orient.getTime(), totalT);
      }
      std::cout << sForm("Total", timer_tot.getTime(), totalT)
		<< std::setw(70) << std::setfill('-') << '-' << std::endl
		<< std::setfill(' ');
    }

    //
    // DEBUG
    //
    if (VERBOSE>4) {
      std::vector<int> levpop(multistep+1, 0), levtot(multistep+1, 0);
      for (auto c : comp->components) {
#if HAVE_LIBCUDA==1
	if (use_cuda) {
	  std::vector<int> clev = c->get_level_lists_cuda();
	  for (int n=0; n<=multistep; n++) levpop[n] += clev[n];
	} else {
	  for (int n=0; n<=multistep; n++) levpop[n] += c->levlist[n].size();
	}
#else
	for (int n=0; n<=multistep; n++) levpop[n] += c->levlist[n].size();
#endif
      }
      
      MPI_Reduce(&levpop[0], &levtot[0], multistep+1, MPI_INT, MPI_SUM, 0,
		 MPI_COMM_WORLD);

      if (myid==0) {
	std::cout << std::endl
		  << std::setw(70) << std::setfill('-') << '-' << std::endl
		  << std::setw(70) << left << "--- Level info" << std::endl
		  << std::setw(70) << std::setfill('-') << '-' << std::endl 
		  << std::setfill(' ') << right;
	std::cout << std::setw(4) << std::left << myid << ": ";
	for (int m=0; m<=multistep; m++) std::cout << std::setw(10) << levpop[m];
	std::cout << std::endl;
	for (int n=1; n<numprocs; n++) {
	  MPI_Recv(levpop.data(), levpop.size(), MPI_INT, n, 191,
		   MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  std::cout << std::setw(4) << std::left << n << ": ";
	  for (int m=0; m<=multistep; m++) std::cout << std::setw(10) << levpop[m];
	  std::cout << std::endl;
	}
      } else {
	MPI_Send(levpop.data(), levpop.size(), MPI_INT, 0, 191, MPI_COMM_WORLD);
      }

      if (myid==0) {
	std::cout << std::endl;
	std::cout << std::setw(4) << std::left << "T" << ": ";
	for (int m=0; m<=multistep; m++) std::cout << std::setw(10) << levtot[m];
	std::cout << " [" << std::accumulate(levtot.begin(), levtot.end(), 0)
		  << "]" << std::endl
		  << std::setw(70) << std::setfill('-') << '-' << std::endl
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
    timer_out  .reset();
    timer_lev  .reset();
    timer_rpt  .reset();
    timer_bal  .reset();
    timer_tot  .reset();
    if (use_cuda) comp->timer_cuda.reset();
    if (use_cuda) comp->timer_orient.reset();
  }

#ifdef USE_GPTL
  GPTLstop("dostep");
#endif
}
