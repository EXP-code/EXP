#ifndef _OutPSN_H
#define _OutPSN_H

#include <OutPSN.H>

/** Write phase-space dumps at regular intervals.
    Each %dump is written into a new file label as <code>filename.n</code>
    where is begins at <code>nbeg</code> and incremented by 1 after each
    file is written.

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
    @param nbeg is suffix of the first phase space %dump
    @param timer set to true turns on wall-clock timer for PS output
*/
class OutPSN : public Output
{

private:

  std::string filename;
  bool timer;
  int nbeg;
  void initialize(void);

public:

  //! Constructor
  OutPSN(const YAML::Node & conf);

  //! Provided by derived class to generate some output
  /*!
    \param nstep is the current time step used to decide whether or not
    to %dump
    \param last should be true on final step to force phase space %dump
    indepentently of whether or not the frequency criterion is met
    \param timer set to true turns on wall-clock timer for PS output
  */
  void Run(int nstep, bool last);

};

#endif
