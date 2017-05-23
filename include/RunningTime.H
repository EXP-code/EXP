#ifndef _RUNNINGTIME_H
#define _RUNNINGTIME_H

#include <deque>
#include <Timer.h>

/*
  A stopwatch that averages over the last N trials
*/
class RunningTime
{
private:
  // The time trial history array
  std::deque<double> boxcar;

  // The current average, the last time, the time scale prefactor
  double current, value, tfac, last, accum;

  // Previous number of states to keep for average
  unsigned nkeep;

  Timer timer;
  bool micro;

  // Add the current value to the history
  void append()
  {
    last   = tfac*timer.getTime().getTotalTime();
    accum += last;
    value  = accum;
    accum  = 0.0;
    boxcar.push_back(value);
    if (boxcar.size() > nkeep) boxcar.pop_front();

    // Compute the average
    current = 0.0;
    for (deque<double>::iterator it=boxcar.begin(); it!=boxcar.end(); it++)
      current += *it;
    current /= boxcar.size();
  }

  // Add to accumulated time
  void lap()
  {
    last   = tfac*timer.getTime().getTotalTime();
    accum += last;
  }

  // Initialize the stop watch
  void initialize(unsigned n, bool microseconds) 
  {
    nkeep = max<unsigned>(n, 1);
    micro = microseconds;
    current = 0.0;
    value   = 0.0;
    accum   = 0.0;
    if (micro) {
      timer.Microseconds();
      tfac = 1.0e-6;
    } else {
      timer.Seconds();
      tfac = 1.0;
    }
  }

public:

  RunningTime()                            { initialize(30, true); }
  RunningTime(unsigned n, bool micro=true) { initialize(n, micro); }

  void Start()
  {
    timer.start();
  }

  void Stop(bool final=true)
  {
    timer.stop();
    if (final) append();
    else       lap();
    timer.reset();
  }

  double getTavg() { return current; }
  double getTime() { return value;   }
  double getLast() { return last;    }
};

#endif
