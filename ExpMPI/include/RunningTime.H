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
  double current, value, tfac;

  // Previous number of states to keep for average
  unsigned nkeep;

  Timer timer;
  bool micro;

  // Add the current value to the history
  void append()
  {
    value = tfac*timer.getTime().getTotalTime();
    boxcar.push_back(value);
    if (boxcar.size() > nkeep) boxcar.pop_front();

    // Compute the average
    current = 0.0;
    for (deque<double>::iterator it=boxcar.begin(); it!=boxcar.end(); it++)
      current += *it;
    current /= boxcar.size();
  }

  // Initialize the stop watch
  void initialize(unsigned n, bool microseconds) 
  {
    nkeep = max<unsigned>(n, 1);
    micro = microseconds;
    current = 0.0;
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

  void Start() { timer.start(); }
  void Stop()  { timer.stop(); append(); timer.reset(); }

  double getTime() { return current; }
};

#endif
