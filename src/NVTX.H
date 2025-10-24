#ifndef NVTX_H
#define NVTX_H

#include <memory>

#include "config_exp.h"

#if HAVE_LIBCUDA==1

#if __CUDACC_VER_MAJOR__ < 12
#include <nvToolsExt.h>
#else
#include <nvtx3/nvToolsExt.h>
#endif 

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

typedef std::shared_ptr<nvTracer> nvTracerPtr;

#endif
