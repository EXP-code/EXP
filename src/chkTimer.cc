#include <cstdlib>
#include <cstdio>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

#include <unistd.h>		// for sleep()

#include "expand.H"
#include "chkTimer.H"

#include "config_exp.h"

extern "C" {
  long rem_time(int);
}
				// 10 minute warning
time_t CheckpointTimer::delta = 600;
				// 2 hours between diagnostic messages
time_t CheckpointTimer::diag  = 7200;

CheckpointTimer::CheckpointTimer()
{
  // Initialization

  initial = final = 0;
  last = current = 0;
  start_time = last_diag = time(0);
  mean = var = 0;
  nsteps = 0;
  firstime = true;
}

void  CheckpointTimer::mark()
{
  if (runtime<0.0) return;

  if (firstime && myid==0) {

    // Get and store the current time
    //
    initial = time(0);

    // Look for job manager
    double runtime0;

    try {
      runtime0 = time_remaining();
    } 
    catch (std::string& msg) {
      std::cout << "CheckpointTimer: no job manager, "
		<< msg << std::endl;
      std::cout << "CheckpointTimer: continuing using default runtime=" 
		<< runtime << std::endl;

      runtime0 = runtime;
    }

    // Compute the final time based on input variable
    final   = initial + static_cast<time_t>(floor(runtime0*3600 + 0.5));

    std::string s_initial = ctime(&initial);
    std::string s_final   = ctime(&final);
    std::cout << "----------------------------------------------------"
	      << "------------------" << std::endl
	      << "CheckpointTimer():" << std::endl
	      << "    current time="  << s_initial
	      << "      final time="  << s_final
	      << "----------------------------------------------------"
	      << "------------------" << std::endl;
    firstime = false;
  }
  //
  // END DEBUG
  //
			      
  last    = current;
  current = time(0);
}

std::ostream& operator<<(std::ostream& out, Time const& T)
{
  time_t hr  = T.t/3600;
  time_t min = (T.t - 3600*hr)/60;
  time_t sec = T.t - 3600*hr - 60*min;
    
  return out << setw(3) << hr  << "h " 
	     << setw(3) << min << "m " 
	     << setw(3) << sec << "s ";
}

bool CheckpointTimer::done()
{
  if (runtime<0.0) return false;

  char flg = 0;
  time_t tr = current - last;

  //
  // The root node decides
  //
  if (myid==0) {
    if (last == 0)                           flg = 0;
    else if (time(0) + tr + delta > final)   flg = 1;
    // Zero otherwise
  }

  MPI_Bcast(&flg, 1, MPI_CHAR, 0, MPI_COMM_WORLD);


  //
  // Diagnostics
  //
  if (myid==0 && last>0) {

    if (nsteps == 0) {
      nsteps++;			// First datum (nsteps=1)
      var = 0.0;
      mean = tr;
    } else {
      nsteps++;			// Update variance and mean (nsteps>1)

      var = ( var*(nsteps-2) + (tr - mean)*(tr - mean)*(nsteps-1)/nsteps ) /
	(nsteps - 1);
      mean = (mean*(nsteps-1) + tr)/nsteps;
    }

    time_t right_now = time(0);
    if (right_now - last_diag > diag) {

      cout << endl
	   << "-------------------------------------------"
	   << "---------------------------"        << endl
	   << "--- Checkpoint timer info -----------------"
	   << "---------------------------"        << endl
	   << "-------------------------------------------"
	   << "---------------------------"        << endl
	   << "Last step: " << Time(tr)            << endl
	   << "Remaining: " << Time(final-time(0)) << endl
	   << "Mean step: " << Time(mean)          << endl
	   << "Root var:  " << Time(sqrt(var))     << endl
	   << "-------------------------------------------"
	   << "---------------------------"        << endl;
      last_diag = right_now;
    }
  }

  return flg ? true : false;
}

std::string CheckpointTimer::exec(std::string& cmd) 
{
  FILE* pipe = popen(cmd.c_str(), "r");
  if (!pipe) return "ERROR";

  char buffer[128];
  std::string result;

  // Give it some time . . .
  sleep(3);

  while(!feof(pipe)) {
    if(fgets(buffer, 128, pipe) != NULL)
      // DEBUG
      if (1) cout << "CheckpointTimer: buffer=" << buffer << endl;
      result += buffer;
  }
  pclose(pipe);

  return result;
}


double CheckpointTimer::time_remaining()
{
  double ret = 0.0;

#ifdef HAVE_LIBSLURM

  std::string env_slurm("SLURM_JOB_ID");
  if (getenv(env_slurm.c_str()) == 0) {
    std::ostringstream sout;
    sout << "No environment variable: SLURM_JOB_ID. Node=" << myid;
    throw GenericError(sout.str(), __FILE__, __LINE__, 1003, true);
  }
  
  std::cout << "----------------------------------------------------"
	    << "------------------" << std::endl
	    << "CheckpointTimer():" << std::endl
	    << "    SLURM detected" << std::endl
	    << "----------------------------------------------------"
	    << "------------------" << std::endl;
  
  std::string job = getenv(env_slurm.c_str());
  long        rem = rem_time(atoi(job.c_str()));

  if (rem<0) {
    std::ostringstream sout;
    sout << "Error obtaining job id from Slurm for node=" << myid;
    throw GenericError(sout.str(), __FILE__, __LINE__, 1003, true);
  }

  ret = static_cast<double>(rem)/3600.0;

#else
  // Specified runtime in hours minus time so far in hours
  ret = runtime - (time(0) - start_time)/3600.0;
#endif

  return ret;
}
