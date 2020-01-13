#ifndef _OutputContainer_H
#define _OutputContainer_H

#include <list>
#include <Output.H>

//! Container class to hold the output chain of container methods 
class OutputContainer
{
private:

  //! Last run time
  double last;

public:
  //! List of output methods currently active
  list<Output*> out;

  //! Constructor
  OutputContainer();

  //! Destructor
  ~OutputContainer();

  //! Will read configuration stanza and initialize requested output methods
  void initialize();

  //! Execute the all methods in the container
  void Run(int n, bool last=false);
};

#endif

