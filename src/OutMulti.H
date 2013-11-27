#ifndef _OutMulti_H
#define _OutMulti_H

#include <OutMulti.H>

/** Dump histogram of multistep levels to file in GNUPLOT format

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
*/
class OutMulti : public Output
{

private:

  string filename;
  void initialize(void);

public:

  //! Constructor
  OutMulti(string& line);

  //! Provided by derived class to generate some output
  /*!
    \param nstep is the current time step used to decide whether or not
    to dump
    \param last should be true on final step to force phase space dump
    indepentently of whether or not the frequency criterion is met
  */
  void Run(int nstep, bool last);

};

#endif
