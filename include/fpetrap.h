/** Some FPE trapping stuff for debugging
 */

#ifndef FPE_TRAP_H
#define FPE_TRAP_H

#include <iostream>
#include <iomanip>
#include <cfenv>
#include <list>

#include <signal.h>
#ifdef HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif

//! Very lean error handler.  Only good for setting a debugger
//! breakpoint.
void my_fpu_handler(int err);

//! Turns on exceptions for invalid, div by zero and overflow, other
//!  bits default.  This will only work for x86 architecture.
void set_fpu_handler(void);

#endif
