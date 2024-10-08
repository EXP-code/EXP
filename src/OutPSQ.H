#ifndef _OutPSQ_H
#define _OutPSQ_H

#include <OutPSQ.H>

/** Write phase-space dumps at regular intervals from each node 
    in component pieces.  These pieces may be reassembled from the info
    in an "instruction" file.

    The assembly instructions for each %dump is written into a new
    file label as <code>filename.n</code> where is begins at
    <code>nbeg</code> and incremented by 1 after each file is written.

    The component pieces are indexed from 0 onward.

    Sending the root process a SIGHUP will cause the first of OutPS,
    OutPSP, OutPSN, OutPSQ, OutCHKPT, or OutCHKPTQ in the Output list
    to execute.  This may not always be possible to signal for batch
    scheduled jobs.

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
    @param nbeg is suffix of the first phase space %dump
    @param timer set to true turns on wall-clock timer for PS output
    @param threads number of threads for binary writes

*/
class OutPSQ : public Output
{

private:

  std::string filename;
  bool real4, timer;
  int nbeg, threads;
  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  OutPSQ(const YAML::Node & conf);

  //! Provided by derived class to generate some output
  /*!
    \param nstep is the current time step used to decide whether or not
    to %dump
    \param mstep is the current multistep level to decide whether or not to dump multisteps
    \param last should be true on final step to force phase space %dump
    indepentently of whether or not the frequency criterion is met
    \param timer set to true turns on wall-clock timer for PS output
    \param threads is the thread count for binary writes
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
