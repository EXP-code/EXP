#ifndef _OutPSR_H
#define _OutPSR_H

#include <OutPSR.H>

/** Write phase-space dumps at regular intervals from each node 
    in component pieces.  These pieces may be reassembled from the info
    in an "instruction" file.

    The assembly instructions for each %dump is written into a new
    file label as <code>filename.n</code> where is begins at
    <code>nbeg</code> and incremented by 1 after each file is written.

    The component pieces are indexed from 0 onward.

    Sending the root process a SIGHUP will cause the first of OutPS,
    OutPSP, OutPSN, OutPSQ, OutPSR, OutCHKPT, or OutCHKPTQ in the
    Output list to execute.  This may not always be possible to signal
    for batch scheduled jobs.

    This class is variant of the OutPSQ class.  This uses OpenMP to
    gang multiple write threads per process.  Whether this increases
    the full bandwidth depends on the details of the file system.
    Testing at Simons suggests that this routine can improve
    throughput for for GPU runs for a small number of processes per
    node but not for many proesses per node.

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
    @param nbeg is suffix of the first phase space %dump
    @param timer set to true turns on wall-clock timer for PS output
    @param threads are the number of writing threads
*/
class OutPSR : public Output
{

private:

  std::string filename;
  bool real4, timer;
  int nbeg, nth;
  void initialize(void);

public:

  //! Constructor
  OutPSR(const YAML::Node & conf);

  //! Provided by derived class to generate some output
  /*!
    \param nstep is the current time step used to decide whether or not
    to %dump
    \param mstep is the current multistep level to decide whether or not to dump multisteps
    \param last should be true on final step to force phase space %dump
    indepentently of whether or not the frequency criterion is met
    \param timer set to true turns on wall-clock timer for PS output
    \param threads are the number of writing threads
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
