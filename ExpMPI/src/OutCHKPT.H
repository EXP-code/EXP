#ifndef _OutCHKPT_H
#define _OutCHKPT_H

#include <OutCHKPT.H>

/** Writes a checkpoint file at regular intervals

    Previously written checkpoint file will be rename to
    <filename>.bak

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
*/
class OutCHKPT : public Output
{

private:

  string filename;
  void initialize(void);

public:

  //! Constructor
  OutCHKPT(string& line);

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
