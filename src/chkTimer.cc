#include <cstdlib>
#include <cstdio>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

#include <unistd.h>		// for sleep()

#include <expand.H>
#include <chkTimer.H>

#include <config.h>

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
  last_diag = time(0);
  mean = var = 0;
  nsteps = 0;
  firstime = true;
}

void  CheckpointTimer::mark()
{
  if (runtime<0.0) return;

  //
  // DEBUG
  //
  if (firstime && myid==0) {

    // Get and store the current time
    initial = time(0);

    

    // Look for job manager
    double runtime0;

    try {
      runtime0 = time_remaining();
    } 
    catch (string& msg) {
      cout << "CheckpointTimer: problem finding job manager, " << msg << endl;
      cout << "CheckpointTimer: continuing using default runtime=" 
	   << runtime << endl;

      runtime0 = runtime;
    }

    // Compute the final time based on input variable
    final   = initial + static_cast<time_t>(floor(runtime0*3600 + 0.5));

    string s_initial = ctime(&initial);
    string s_final   = ctime(&final);
    cout << "----------------------------------------------------"
	 << "------------------" << endl
	 << "CheckpointTimer():" << endl
	 << "    current time="  << s_initial
	 << "      final time="  << s_final
	 << "----------------------------------------------------"
	 << "------------------" << endl;
    firstime = false;
  }
  //
  // END DEBUG
  //
			      
  last    = current;
  current = time(0);
}

ostream& operator<<(ostream& out, Time const& T)
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

string CheckpointTimer::exec(string& cmd) 
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

  string env_slurm("SLURM_JOB_ID");
  if (getenv(env_slurm.c_str()) == 0)
    throw string("No environment variable: SLURM_JOB_ID");
  
  std::cout << "----------------------------------------------------"
	    << "------------------" << endl
	    << "CheckpointTimer():" << endl
	    << "    SLURM detected" << endl
	    << "----------------------------------------------------"
	    << "------------------" << endl;

  string job = getenv(env_slurm.c_str());
  long   rem = rem_time(atoi(job.c_str()));

  if (rem<0)
    throw string("Error obtaining job id from slurm");

  ret = static_cast<double>(rem)/3600.0;

  return ret;

#else

  string env_pbs("PBS_JOBID");

  if (getenv(env_pbs.c_str()) == 0)
    throw string("No environment variable: PBS_JOBID");

  std::cout << "----------------------------------------------------"
	    << "------------------" << endl
	    << "CheckpointTimer():" << endl
	    << "    PBS detected"   << endl
	    << "----------------------------------------------------"
	    << "------------------" << endl;

  string job = getenv(env_pbs.c_str());

  string command = "showstart " + job;
  string result = exec(command);

  istringstream sin(result);
  char line[128];
  string tag("Earliest completion in ");
  bool found = false;

  while (sin) {
    sin.getline(line, 128);
    string sline(line);
    
    if (sline.find(tag) != string::npos) {

      found = true;

      string chop = sline.substr(sline.find(tag)+tag.length());
      chop = chop.substr(0, chop.find("on"));
      chop = chop.substr(chop.find_first_not_of(" "));
      chop = chop.substr(0, chop.find(" "));

#ifdef DEBUG
      cout << "String: " << chop << endl;
#endif

      vector<string> items;
      int days=0, hours=0, mins=0, secs=0;
      
      while (chop.size()) {
	if (chop.find(":") != string::npos) {
	  items.push_back(chop.substr(0, chop.find(":")));
	} else {
	  items.push_back(chop);
	  break;
	}
	chop = chop.substr(chop.find(":")+1);
      }

#ifdef DEBUG
      for (unsigned n=0; n<items.size(); n++)
	cout << setw(3) << n << " <" << items[n] << ">" << endl;
#endif

      if (items.size()) {
	secs = atoi(items.back().c_str());
	items.pop_back();
      }

      if (items.size()) {
	mins = atoi(items.back().c_str());
	items.pop_back();
      }

      if (items.size()) {
	hours = atoi(items.back().c_str());
	items.pop_back();
      }

      if (items.size()) {
	days = atoi(items.back().c_str());
	items.pop_back();
      }

      ret = static_cast<double>(((days*24 + hours)*60 + mins)*60 + secs) / 3600.0;

      break;
    }

  }

  if (!found) 
    throw string("PBS command failure, no PBS?? ") +
      "[" + command + "]==>[" + result + "]";

#endif

  return ret;
}
