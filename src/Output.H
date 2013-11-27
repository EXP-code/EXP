#ifndef _Output_H
#define _Output_H

#include <string>
#include <vector>
#include <utility>
#include <StringTok.H>

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
*/
class Output
{

private:

protected:

  //! Parameter database
  list< pair<string, string> > namevalue;

  //! Output frequency (should be used by derived classes)
  int nint;

  //! Used to assign additional specific parameters, create storage, etc.
  virtual void initialize(void) = 0;

public:

  //! Id string
  string id;

  //! Constructor
  Output(string& line);

  //! Destructor
  virtual ~Output() {}

  //! Provided by derived class to generate some output
  virtual void Run(int nstep, bool final) = 0;

  //! Get token by name from initialization string
  int get_value(const string& name, string& value);

  //! Returns a map from token(n)=value pairs
  std::map<int, string> get_value_array(const string& name);

  //! Error
  void bomb(const string& msg);
  
};

#endif
