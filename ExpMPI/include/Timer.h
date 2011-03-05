// This is really -*- C++ -*-

#ifndef _TIMER_H
#define _TIMER_H

#include <iostream>
#include <cstdlib>

#include <time.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>


//! This class is a container for returning Timer results
class TimeElapsed
{
private:
  long userTime, systemTime, realTime;
  bool micro;

public:

  /** Setting this to true returns the wallclock time in the () operator
   for TimeElapsed rather than the system + user time for the process.
   Obviously, this is only for diagnostic and unlikely to affect an
   simulation values (default value: false). */
  static bool use_real;

  //! Null constructor
  TimeElapsed() : userTime(0), systemTime(0), realTime(0), micro(false) {}

  //! Copy constructor
  TimeElapsed(const TimeElapsed& p) {
    userTime   = p.userTime;
    systemTime = p.systemTime;
    realTime   = p.realTime;
    micro      = p.micro;
  }

  //! Standard constructor
  TimeElapsed(long user, long system, long real, bool type) {
    userTime   = user;
    systemTime = system;
    realTime   = real;
    micro      = type;
  }

  //! True if raw long values are in seconds
  bool usingSeconds() { return !micro; }

  //! True if raw long values are in microseconds
  bool usingMicroseconds() { return micro; }

  //! Return user time (in seconds or microseconds)
  long getUserTime() {
    return userTime;
  }
  
  //! Return system time (in seconds or microseconds)
  long getSystemTime() {
    return systemTime;
  }
  
  //! Return total tiem (in seconds or microseconds)
  long getTotalTime() {
    return userTime + systemTime;
  }
      
  //! Get the wallclock time (in seconds or microseconds)
  long getRealTime() {
    return realTime;
  }
  
  /** Return the time in seconds always.  The precision is determined
      by the Timer instance and the type of return value is set by the
      use_real static variable */
  double operator()() { 
    double ret = 0.0;
    if (use_real) ret = realTime;
    else          ret = userTime + systemTime;
    if (micro)    ret *= 1.0e-6;
    return ret;
  }
};


//! This is a diagnostic timer (wraps up the POSIX time calls)
class Timer
{
private:

  timeval begin, end; 
  struct rusage beginusage, endusage;

  //@{
  //! the time counted so far
  long utime, stime, rtime;
  //@}
  
  //! indicate whether the timer has been started
  bool started;

  /** if true then return number of microseconds for user and system
      time; else return number of seconds */
  bool precision; 

 public:
  /** Main constructor: if precision is set to true, then the real and
      system time will be measured in microseconds */
  Timer(bool precision = false) {
    utime = stime = rtime = 0;
    started = false;
    this->precision = precision;
  }

  //! Copy consructor
  Timer(const Timer& t) {
    begin      = t.begin;
    end        = t.end;
    beginusage = t.beginusage;
    endusage   = t.endusage;
    utime      = t.utime;
    stime      = t.stime;
    rtime      = t.rtime;
    started    = t.started;
    precision  = t.precision;
  }

  //! Set precision to seconds
  void Seconds()      { precision = false; }

  //! Set precision to seconds
  void Microseconds() { precision = true; }

  //! Return precision (true = microseconds, false = seconds)
  bool Precision()    { return precision; }

  //! start timer, if already started then do nothing
  void start() {
    if (started) return;
    started = true;

    if (gettimeofday(&begin, NULL))
      std::cerr << "gettimeofday error!";

    if (getrusage(RUSAGE_SELF, &beginusage) == -1)
      std::cerr << "getrusage error!";
  }
  
  /** stop timer and return time measured so far.  if timer was
      stopped the time returned will be 0. */
  TimeElapsed stop() {
    if (!started)
      return TimeElapsed(0, 0, 0, precision);
    
    if (gettimeofday(&end, NULL))
      std::cerr << "gettimeofday error!";

    if (getrusage(RUSAGE_SELF, &endusage) == -1)
      std::cerr << "getrusage error!";

    started = false;

    if (precision) {
      long uusecs = 
	(endusage.ru_utime.tv_sec   - beginusage.ru_utime.tv_sec) * 1000000 +
	(endusage.ru_utime.tv_usec  - beginusage.ru_utime.tv_usec);

      utime += uusecs;
      
      long susecs = 
	(endusage.ru_stime.tv_sec   - beginusage.ru_stime.tv_sec) * 1000000 +
	(endusage.ru_stime.tv_usec  - beginusage.ru_stime.tv_usec);

      stime += susecs;
      
      long rusecs = 
	(end.tv_sec - begin.tv_sec) * 1000000 + end.tv_usec - begin.tv_usec;

      rtime += rusecs;
      
    } else {
      long usecs = (endusage.ru_utime.tv_sec - beginusage.ru_utime.tv_sec);

      utime += usecs;
	  
      long ssecs = (endusage.ru_stime.tv_sec - beginusage.ru_stime.tv_sec);

      stime += ssecs;
      
      long rsecs = 
	(end.tv_sec  - begin.tv_sec ) + 
	(end.tv_usec - begin.tv_usec)/1000000;

      rtime += rsecs;
    }
    
    return TimeElapsed(utime, stime, rtime, precision);
  }

  /** reset the timer, this will reset the time measured to 0 and will
      leave the timer in the same status (started or stopped). */
  void reset() {
    utime = stime = rtime = 0;

    if (started) {
      if (gettimeofday(&begin, NULL))
	std::cerr << "gettimeofday error!";
      
      if (getrusage(RUSAGE_SELF, &beginusage) == -1)
	std::cerr << "getrusage error!";
    }
  }

  //! return time measured up to this point.
  TimeElapsed getTime() {
    if (!started)
      return TimeElapsed(utime, stime, rtime, precision);
    
    if (gettimeofday(&end, NULL))
      std::cerr << "gettimeofday error!";
    
    if (getrusage(RUSAGE_SELF, &endusage) == -1)
      std::cerr << "getrusage error!";

    if (precision) {
      long uusecs = 
	(endusage.ru_utime.tv_sec  - beginusage.ru_utime.tv_sec ) * 1000000 +
	(endusage.ru_utime.tv_usec - beginusage.ru_utime.tv_usec);
      
      long susecs = 
	(endusage.ru_stime.tv_sec  - beginusage.ru_stime.tv_sec ) * 1000000 +
	(endusage.ru_stime.tv_usec - beginusage.ru_stime.tv_usec);
      
      long rusecs = 
	(end.tv_sec - begin.tv_sec) * 1000000 + end.tv_usec - begin.tv_usec;
      
      return TimeElapsed(utime + uusecs, 
			 stime + susecs, 
			 rtime + rusecs,
			 precision);
    } else {
      long usecs = (endusage.ru_utime.tv_sec - beginusage.ru_utime.tv_sec);
      
      long ssecs = (endusage.ru_stime.tv_sec - beginusage.ru_stime.tv_sec);
      
      long rsecs = 
	(end.tv_sec  - begin.tv_sec ) +
	(end.tv_usec - begin.tv_usec)/1000000;
      
      return TimeElapsed(utime + usecs, 
			 stime + ssecs, 
			 rtime + rsecs,
			 precision);
    }
  }
  
  //! Return the status of the timer
  bool isStarted() { return started; }

};

#endif // _TIMER_H
