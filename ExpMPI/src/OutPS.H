#ifndef _OutPS_H
#define _OutPS_H

#include <OutPS.H>

/** Append phase space dumps to specified file at regular intervals

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
*/
class OutPS : public Output
{

private:

  string filename;
  void initialize(void);

public:

  //! Constructor
  OutPS(string& line);

  //! Provided by derived class to generate some output
  /*!
    <param> nstep is the current time step used to decide whether or not
    to dump
    <param> last should be true on final step to force phase space dump
    indepentently of whether or not the frequency criterion is met
  */
  void Run(int nstep, bool last);

};

#endif
