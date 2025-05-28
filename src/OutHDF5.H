#ifndef _OutHDF5_H
#define _OutHDF5_H

#include <OutHDF5.H>

/** Write phase-space dumps at regular intervals from each node in
    component pieces.  This writer is indended to be Gadget-4 HDF5
    compliant.

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
    @param nbeg is suffix of the first phase space %dump
    @param timer set to true turns on wall-clock timer for PS output
    @param threads number of threads for binary writes
    @param noids set to true turns off particle id writing

*/
class OutHDF5 : public Output
{

private:

  std::string filename;
  bool real4=true, real8=false, ids=true, timer;
  int nbeg, threads;
  void initialize(void);
  std::vector<double> masses;
  std::vector<bool> multim;

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

  //! Check for single particle mass on the first invocation
  void checkParticleMasses();

public:

  //! Constructor
  OutHDF5(const YAML::Node & conf);

  //! Provided by derived class to generate some output
  /*!
    \param nstep is the current time step used to decide whether or not
    to %dump
    \param mstep is the current multistep level to decide whether or not to dump multisteps
    \param last should be true on final step to force phase space %dump
    indepentently of whether or not the frequency criterion is met
    \param timer set to true turns on wall-clock timer for PS output
    \param threads is the thread count for binary writes
    \param noids set to true turns off particle id writing
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
