#ifndef _Centering_H_
#define _Centering_H_

#include <ParticleReader.H>

namespace Utility
{
  std::vector<double> getDensityCenter(PR::PRptr reader, int Ndens=32, int stride=1);
  std::vector<double> getCenterOfMass (PR::PRptr reader);
}

#endif
