
/** Some FPE trapping stuff for debugging
 */

#ifndef FPE_TRAP_H
#define FPE_TRAP_H

#include <iostream>
#include <iomanip>
#include <list>

#include <fenv.h>
#include <signal.h>
#ifdef HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif

/// Global to hold error code
volatile int my_err;

/// Very lean error handler.  Only good for setting a debugger breakpoint.
static void my_fpu_handler(int err)
{
  my_err = err;
  std::cerr << "FP trapped error=" << my_err << std::endl;
}

static void (*oldhandler)(int);	// keep a reference to the initial value
				// so it can be restored later

/**
   Turns on exceptions for invalid, div by zero and overflow, other  
   bits default.  This will only work for x86 architecture.
 */
void set_fpu_handler(void)
{
#ifdef HAVE_FPU_CONTROL_H 
				// Set control flag (see fpu_control.h)
  short cw = 0x0372;
  _FPU_SETCW(cw);
#endif
  
  oldhandler = signal(SIGFPE, my_fpu_handler);
  if (SIG_ERR == oldhandler) {
    std::cerr << "EXP: Cannot install floating-point exception handler" << std::endl;
    exit(-1);
  } else {
    std::cout << "EXP: Floating point-exception handler installed" << std::endl;
  }

#ifdef HAVE_FPU_CONTROL_H
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
  //
  // Print enabled flags to root node
  //
  if (myid==0) {
    const std::list<std::pair<int, std::string>> flags =
      {	{FE_DIVBYZERO, "divide-by-zero"},
	{FE_INEXACT,   "inexact"},
	{FE_INVALID,   "invalid"},
	{FE_OVERFLOW,  "overflow"},
	{FE_UNDERFLOW, "underflow"} };
    
    int _flags = fegetexcept();
    std::cout << "---- Enabled FE flags: <";
    for (auto v : flags) {
      if (v.first & _flags) std::cout << v.second << ' ';
    }
    std::cout << "\b>" << std::endl;
  }
#endif
  
}

#endif
