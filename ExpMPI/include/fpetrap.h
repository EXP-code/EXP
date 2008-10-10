
/** Some FPE trapping stuff for debugging
 */

#ifndef FPU_TRAP_H
#define FPE_TRAP_H

#include <signal.h>
#include <fpu_control.h>

/// Global to hold error code
volatile int my_err;

/// Very lean error handler.  Only good for setting a debugger breakpoint.
static void my_fpu_handler(int err)
{
  my_err = err;
  cerr << "FP trapped error=" << my_err << endl;
}

static void (*oldhandler)(int);	// keep a reference to the initial value
				// so it can be restored later

/**
   Turns on exceptions for invalid, div by zero and overflow, other  
   bits default.  This will only work for x86 architecture.
 */
void set_fpu_handler(void)
{
				// Set control flag (see fpu_control.h)
  short cw = 0x0372;
  _FPU_SETCW(cw);
  
  oldhandler = signal(SIGFPE, my_fpu_handler);
  if (SIG_ERR == oldhandler) {
    cerr << "cannot install floating point exception handler";
    exit(-1);
  }
}

#endif
