// Example of using the version string and parsed octet of integers from libvars.H

#include <iostream>
#include "libvars.H"
using namespace __EXP__;

int main()
{
  std::cout << "Version string: " << VERSION << '\n';
  std::cout << "Parsed octet of integers:\n";
  std::cout << "-- Major=" << exp_build.major << '\n';
  std::cout << "-- Minor=" << exp_build.minor << '\n';
  std::cout << "-- Patch=" << exp_build.patch << '\n';
  return 0;
}
