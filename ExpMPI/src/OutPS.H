#ifndef _OutPS_H
#define _OutPS_H

#include <OutPS.H>

//! Print field quantities for spherical and cylindrical bases
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
