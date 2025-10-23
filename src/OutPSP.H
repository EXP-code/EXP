#ifndef _OutPSP_H
#define _OutPSP_H

#include "OutPSP.H"

/** Write phase-space dumps at regular intervals using MPI-IO

    Each %dump is written into a new file label as
    <code>filename.n</code> where is begins at <code>nbeg</code> and
    incremented by 1 after each file is written.  This implementation
    bunches particles into memory and uses a collective MPI_File write
    to push the bunch to disk at the appropriate offset.  The MPI-IO
    implementation does not maintain the particle sequence so
    indexing=1 needs to be added as a Component option if the particle
    index is needed for post processing.  This version is much faster
    than OutPSN and OutCHKPT which uses MPI gather to the root node
    before writing.  This is the newest version for PSP output and has
    not be extensively tested as of 10/2018.

    Sending the root process a SIGHUP will cause the first of OutPS,
    OutPSP, OutPSN, OutPSQ, OutCHKPT, or OutCHKPTQ in the output list
    to execute.  This may not always be possible to signal for batch
    scheduled jobs.

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
    @param nbeg is suffix of the first phase space %dump
    @param real4 indicates floats for real PS quantities
    @param nagg is the number of MPI-IO aggregators
*/
class OutPSP : public Output
{

private:

  std::string filename, nagg;
  bool real4, timer;
  int nbeg;

  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  OutPSP(const YAML::Node& conf);

  //! Provided by derived class to generate some output
  /*!
    \param nstep is the current time step used to decide whether or not
    to %dump
    \param mstep is the current multistep level to decide whether or not to dump multisteps
    \param last should be true on final step to force phase space %dump
    indepentently of whether or not the frequency criterion is met
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
