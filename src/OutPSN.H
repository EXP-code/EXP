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
*/
class OutPSN : public Output
{

private:

  string filename;
  int nbeg;
  void initialize(void);

public:

  //! Constructor
  OutPSN(string& line);

  //! Provided by derived class to generate some output
  /*!
    \param nstep is the current time step used to decide whether or not
    to %dump
    \param last should be true on final step to force phase space %dump
    indepentently of whether or not the frequency criterion is met
  */
  void Run(int nstep, bool last);

};

#endif
