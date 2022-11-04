#ifndef _CheckpointTimer_H
#define _CheckpointTimer_H

#include <time.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//! A "helper" class to print time_t intervals
class Time
{
protected:
  time_t t;
    
public:
  //! Constructor for time_t input
  Time(time_t T) { t = T; }

  //! Constructor for double input
  Time(double T) { t = static_cast<time_t>(floor(T+0.5)); }
  
  //! std output stream manipulator
  friend ostream& operator<<(ostream& out, Time const& T);
};


/**
  This class keeps track of time between time steps and calls for
  a checkpoint and signals an exit if there will not be time to reach 
  the next checkpoint
*/
class CheckpointTimer
{
private:

  time_t initial, final, last, current, last_diag, start_time;
  unsigned nsteps;
  double mean, var;
  bool firstime;

  //! Execute a command and pipe the results to string
  string exec(string& cmd);

  //! Get the remaining execution time
  double time_remaining();

public:
  
  //! Safety time interval       [Default: 10 minute warning]
  static time_t delta;

  //! Diagnostic output interval [Default: 2 hours]
  static time_t diag;

  //! Constructor
  CheckpointTimer();

  //! Mark the current time
  void mark();

  //! Time to quit
  bool done();

};

#endif
