#ifndef _OutPS_H
#define _OutPS_H

#include <OutPS.H>

/** Append phase space dumps to specified file at regular intervals

    Sending the root process a SIGHUP will cause the first of OutPS,
    OutPSP, OutPSN, OutPSQ, OutCHKPT, or OutCHKPTQ in the Output list
    to execute.  This may not always be possible to signal for batch
    scheduled jobs.

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
    @param timer set to true turns on wall-clock timer for PS output
*/
class OutPS : public Output
{

private:

  std::string filename;
  bool timer;
  
  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  OutPS(const YAML::Node& conf);

  //! Provided by derived class to generate some output
  /*!
    \param nstep is the current time step used to decide whether or not
    to dump
    \param mstep is the current multistep level to decide whether or not to dump multisteps
    \param last should be true on final step to force phase space dump
    indepentently of whether or not the frequency criterion is met
    \param timer set to true turns on wall-clock timer for PS output
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
