#include <Timer.h>

// Setting this to true returns the wallclock time in the () operator
// for TimeElapsed rather than the system + user time for the process.
// Obviously, this is only for diagnostic and unlikely to affect an
// simulation values.
//
bool TimeElapsed::use_real = false;
