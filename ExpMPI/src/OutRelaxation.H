#ifndef _OutRelaxation_H
#define _OutRelaxation_H

#include <Output.H>

//! Output routine for generateRelaxation ("external force") module
class OutRelaxation : public Output
{

private:

  string fname;

  int epos;
  void initialize(void);

public:

  //! Constructor
  OutRelaxation(string& line);
  //! Destructor
  virtual ~OutRelaxation() {};

  //! Do some output
  /*!
    <param> nstep is the current time step index
    <param> final should be true on last step to force output
  */
  void Run(int nstep, bool final=false);

};

#endif
