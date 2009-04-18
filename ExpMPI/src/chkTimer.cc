#include <expand.h>
#include <chkTimer.H>

				// 10 minute warning
time_t CheckpointTimer::delta = 600;

CheckpointTimer::CheckpointTimer()
{
  // Get and store the current time
  
  initial = time(0);

  // Compute the final time based on input variable

  final   = initial + static_cast<time_t>(floor(runtime*3600 + 0.5));

  // Initialization

  last = current = 0;
  mean = var = 0;
  nsteps = 0;
  firstime = true;
}

void  CheckpointTimer::mark()
{
  //
  // DEBUG
  //
  if (firstime && myid==0) {
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

  MPI_Bcast(&flg, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  //
  // DEBUG
  //
  if (true) {

    if (myid==0 && last>0) {

      if (nsteps == 0) {
	nsteps++;		// First datum (nsteps=1)
	var = 0.0;
	mean = tr;
      } else {
	nsteps++;		// Update variance and mean (nsteps>1)

	var = ( var*(nsteps-2) + (tr - mean)*(tr - mean)*(nsteps-1)/nsteps ) /
	  (nsteps - 1);
	mean = (mean*(nsteps-1) + tr)/nsteps;
      }

      cout << "-------------------------------------------"
	   << "---------------------------" << endl
	   << "--- Checkpoint timer info -----------------"
	   << "---------------------------" << endl
	   << "-------------------------------------------"
	   << "---------------------------" << endl
	   << "Last step: " << Time(tr)            << endl
	   << "Remaining: " << Time(final-time(0)) << endl
	   << "Mean step: " << Time(mean)          << endl
	   << "Root var:  " << Time(sqrt(var))     << endl
	   << "-------------------------------------------"
	   << "---------------------------" << endl;
    }
  }
  //
  // END DEBUG
  //

  return flg ? true : false;
}
