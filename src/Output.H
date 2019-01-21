#ifndef _Output_H
#define _Output_H

#include <string>
#include <vector>
#include <utility>
#include <StringTok.H>
#include <Timer.h>

//! Output routines.
/*! These are designed to be run in a chain at each step.  Each Output
  class is responsible for determining when to run.  Initialization
  parameters are parsed on construction and available to the
  contructor for initializing control variables.

  The class is also responsible for cycling through the component list
  and operating on each phase space accordingly.  NB: the alternative
  would be to allow each component to define an output method.  However,
  this can easily result in a chaos of specialized output streams and so
  I rejected this approach.

  The stopWatch Timer is implemented in the base class but can be used
  by any derived class to time the output event.  See OutPS for an
  example.
*/
class Output
{

private:

protected:

  //! Configuration
  const YAML::Node& conf;

  //! Output frequency (should be used by derived classes)
  int nint;

  //! Use output timer
  bool timer;

  //! Output timer
  Timer stopWatch;

  //! Used to assign additional specific parameters, create storage, etc.
  virtual void initialize(void) = 0;

public:

  //! Id string
  string id;

  //! Constructor
  Output(const YAML::Node& conf);

  //! Destructor
  virtual ~Output() {}

  //! Provided by derived class to generate some output
  virtual void Run(int nstep, bool final) = 0;

  //! Error
  void bomb(const string& msg);
  
};

#endif
