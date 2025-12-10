#ifndef _OutputContainer_H
#define _OutputContainer_H

#include <thread>
#include <list>

#include "Output.H"

//! Container class to hold the output chain of container methods 
class OutputContainer
{
private:

  //! Last run time
  double last;

public:

  //! A list of concurrent threads for asychronous processing
  std::list<std::shared_ptr<std::thread>> cproc;

  //! List of output methods currently active
  std::list<Output*> out;

  //! Constructor
  OutputContainer();

  //! Destructor
  ~OutputContainer();

  //! Will read configuration stanza and initialize requested output methods
  void initialize();

  //! Execute the all methods in the container
  void Run(int nstep, int mstep=std::numeric_limits<int>::max(), bool final=false);
};

#endif

