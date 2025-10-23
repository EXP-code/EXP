#ifndef _RUNNINGTIME_H
#define _RUNNINGTIME_H

#include <deque>
#include "Timer.H"

/*
  A stopwatch that averages over the last N trials
*/
class RunningTime
{
private:
  // The time trial history array
  std::deque<double> boxcar;

  // The current average, the last time, the time scale prefactor
  double current, value, last, accum;

  // Previous number of states to keep for average
  unsigned nkeep;

  Timer timer;

  // Add the current value to the history
  void append()
  {
    last   = timer.getTime();
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
    last   = timer.getTime();
    accum += last;
  }

  // Initialize the stop watch
  void initialize(unsigned n) 
  {
    nkeep = max<unsigned>(n, 1);
    current = 0.0;
    value   = 0.0;
    accum   = 0.0;
  }

public:

  RunningTime()           { initialize(30); }
  RunningTime(unsigned n) { initialize(n ); }

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
