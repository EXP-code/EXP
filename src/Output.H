#ifndef _Output_H
#define _Output_H

#include <string>
#include <vector>
#include <utility>
#include <set>

#include "StringTok.H"
#include "Timer.H"

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

  //! Multistep output frequency (also used by derived classes)
  int nintsub;

  //! Use output timer
  bool timer;

  //! Output timer
  Timer stopWatch;

  //! Used to assign additional specific parameters, create storage, etc.
  virtual void initialize(void) = 0;

  //! nintsub warning message
  void nintsub_warning(const std::string& name)
  {
    if (myid==0) std::cout << name << ": substep output is disabled because "
			   << "the phase space dump may include particles at "
			   << "later times for low-level particles."
			   << std::endl
			   << name << ": please recompile with -DALLOW_NINTSUB "
			   << "if you want this option anyway."
			   << std::endl;
  }

  //! Current YAML keys to check configuration
  std::set<std::string> current_keys;

public:

  //! Id string
  string id;

  //! Constructor
  Output(const YAML::Node& conf);

  //! Destructor
  virtual ~Output() {}

  //! Provided by derived class to generate some output
  virtual void Run(int nstep, int mstep, bool final) = 0;

  //! Return unmatched parameters
  std::set<std::string> unmatched() { return current_keys; }
};

#endif
