#ifndef _PERMUTATION_H_
#define _PERMUTATION_H_

#include <algorithm>
#include <numeric>
#include <random>

/** Make a permutation index.  Share will MPI nodes if MPI is active.

    This implementation should fall back politely to single processor
    use.
*/
class permutation
{
public:

  /** Initialize and generate the identity permutation

      \par Parameters:

      @param n is the vecctor rank

      @param full set to true will seed every value in the mt19937
      state. This will produce better randomness at higher startup
      cost.
  */
  permutation(unsigned n, bool full=false) : perm(n), g(std::random_device{}())
  {
    if (full) seed19937();	// Hard-core mt19937 initialization
    std::iota(perm.begin(), perm.end(), unsigned(0));
  }

  //! Create and share the random permutation
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

  //! Seed the Mersenne Twister from the random device
  void seed19937()
  {
    // A hardware random # generator, hopefully
    std::random_device rdev;

    // Generate an array that fills the state
    std::seed_seq::result_type data[std::mt19937::state_size];

    // Generate random number for the full state
    std::generate_n(data, std::mt19937::state_size, std::ref(rdev));
    
    // Fill the sequence
    std::seed_seq prng_seed (data, data + std::mt19937::state_size);

    // Finally, seed the mt19937
    g.seed(prng_seed);
  }
};

#endif
