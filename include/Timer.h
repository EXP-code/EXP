// This is really -*- C++ -*-

#ifndef _TIMER_H
#define _TIMER_H

#include <iostream>
#include <cstdlib>
#include <chrono>

//! This is a diagnostic timer (wraps up the POSIX time calls)
class Timer
{
private:

  //@{
  //! Time structures
  std::chrono::high_resolution_clock::time_point begin, end; 
  //@}

  //@{
  //! the time counted so far
  double rtime;
  //@}
  
  //! indicate whether the timer has been started
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

  //! start timer, if already started then do nothing
  void start()
  {
    if (started) return;
    started = true;
    begin   = std::chrono::high_resolution_clock::now();
  }
  
  /** stop timer and return time measured so far.  if timer was
      stopped the time returned will be 0. */
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

  /** reset the timer, this will reset the time measured to 0 and will
      leave the timer in the same status (started or stopped). */
  void reset()
  {
    rtime   = 0;
    started = false;
  }

  //! return time measured up to this point.
  double getTime() { return rtime; }
  
  //! Return the status of the timer
  bool isStarted() { return started; }

};

#endif // _TIMER_H
