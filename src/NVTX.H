#ifndef NVTX_H
#define NVTX_H

#include <boost/shared_ptr.hpp>

#include "config.h"

#if HAVE_LIBCUDA==1

#include <nvToolsExt.h>

class nvTracer {
private:
  static unsigned count;

public:
  nvTracer(const char* name);

  ~nvTracer() {
    nvtxRangePop();
  }
};
#else
class nvTracer {
public:
  nvTracer(const char* name) {}
  ~nvTracer() {}
};
#endif

typedef boost::shared_ptr<nvTracer> nvTracerPtr;

#endif
