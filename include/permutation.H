#ifndef _PERMUTATION_H_
#define _PERMUTATION_H_

#include <algorithm>
#include <numeric>
#include <random>

/** Make a permutation index.  Share will MPI nodes if MPI is active.

    Should fall back politely to single processor use.

    This will initialze the Mersenne Twister from the random device.
    This initialization is not cryptographically sound, I know, but
    that doesn't matter here.
*/
struct permutation
{
  permutation(unsigned n) : perm(n), g(std::random_device{}())
  {
    if (myid==0) std::iota(perm.begin(), perm.end(), unsigned(0));
  }

  void shuffle() {
    if (myid==0) std::shuffle(perm.begin(), perm.end(), g);
    if (numprocs>1)
      MPI_Bcast(perm.data(), perm.size(), MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  }

  //! The permutation operator
  size_t operator[](size_t n) const { return perm[n]; }

private:
  std::vector<size_t> perm;
  std::mt19937 g;
};

#endif
