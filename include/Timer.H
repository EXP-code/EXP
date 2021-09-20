#ifndef _TIMER_H
#define _TIMER_H

#include <iostream>
#include <cstdlib>
#include <chrono>

//! This is a diagnostic stopwatch-like timer implemented using std::chrono
class Timer
{
private:

  //@{
  //! Time structures (will initialize to zero by default)
  std::chrono::high_resolution_clock::time_point begin, end; 
  //@}

  //@{
  //! The time measured so far
  double rtime;
  //@}
  
  //! Indicate whether the timer has been started
  bool started;

 public:

  //! Constructor
  Timer()
  {
    rtime   = 0;
    started = false;
  }

  //! Copy consructor
  Timer(const Timer& t)
  {
    begin      = t.begin;
    end        = t.end;
    rtime      = t.rtime;
    started    = t.started;
  }

  //! Start timer, if already started then do nothing
  void start()
  {
    if (started) return;
    started = true;
    begin   = std::chrono::high_resolution_clock::now();
  }
  
  /** Stop timer and return time measured so far.  If timer was
      stopped the time returned will be the current elapsed time. */
  double stop()
  {
    if (started) {
      end     = std::chrono::high_resolution_clock::now();
      started = false;
      std::chrono::duration<double> duration = end - begin;
      rtime  += duration.count();
    }
    return rtime;
  }

  /** Reset the timer, this will reset the time measured to 0 and
      reset the zero point if the timer is active */
  void reset() {
    if (started) {
      begin   = std::chrono::high_resolution_clock::now();
    }
    rtime = 0.0;
  }

  //! Return time measured up to this point.
  double getTime()
  {
    if (started) {
      end     = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration = end - begin;
      rtime  += duration.count();
      begin   = end;
    }
    return rtime;
  }
  
  //! Return time measured up to this point.
  double operator()() { return getTime(); }
  
  //! Return the status of the timer
  bool isStarted() { return started; }

};

#endif // _TIMER_H
