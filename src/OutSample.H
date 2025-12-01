#ifndef _OutSample_H
#define _OutSample_H

#include "OutSample.H"

/** @file OutSample.H

    @brief Dump coefficients at each interval
    
    @param filename of the coefficient file

    @param name of the component to dump
  
    @note This output class requires a different strategy than most of
    the others.  Rather than maintaining its own step counter, it
    checks the component's subsampling status to determine whether
    subsample data is available.  Then it updates the HDF5 file.  The
    `mstep` and `last` parameters in the `Run`.  The component's
    `PotAccel` instance must support subsampling by overriding and
    implementing `subsampleReady()` as well as the underlying
    `getSubsample()` member to retrieve the data.
*/
class OutSample : public Output
{

private:

  std::string filename;
  double prev = -std::numeric_limits<double>::max();
  Component *tcomp;

  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  OutSample(const YAML::Node& conf);

  //! Generate the output
  /*!
    \param nstep is the current time step used to decide whether or not
    to dump
    \param last should be true on final step to force phase space dump
    indepentently of whether or not the frequency criterion is met
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
