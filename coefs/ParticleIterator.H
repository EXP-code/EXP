#ifndef _ParticleIterator_H_
#define _ParticleIterator_H_

#include <ParticleReader.H>

namespace Utility
{
  //! Callback function signature for user computation over phase
  //! space
  using Callback =
    std::function<void(
		       double,		     // mass
		       std::vector<double>&, // position
		       std::vector<double>&, // velocity
		       uint64_t		     // index
		       )>;
  
  //! Apply the callback to all particles
  void particleIterator(PR::PRptr reader, const Callback& func);
}

#endif
