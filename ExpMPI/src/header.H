#ifndef header_H
#define header_H

#include <iostream>

using namespace std;

//! Describe the phase space dump
struct MasterHeader {
  //! Current time for phase space
  double time;
  //! Number of particles in entire phase space
  int ntot;
  //! Number of individual components
  int ncomp;
};

//! Describe a component in a phase space dump
class ComponentHeader 
{
public:
  //! Number of bodies in the component
  int nbod;

  //! Number of integer attributes per particle
  int niatr;

  //! Number of real (double) attributes per particle
  int ndatr;

  //! Size of info string
  int ninfochar;

  //! Info string
  char *info;
  
  //! Default size of info string (if not specified explicitly or read)
  static int defaultInfoSize;

  //! Default constructor: uses defaultInfoSize by default
  ComponentHeader() 
  {
    ninfochar = defaultInfoSize;
    info = new char [ninfochar];
  }

  //! Construct component with info field of length n
  ComponentHeader(int n)
  {
    ninfochar = n;
    info = new char [ninfochar];
  }

  //! Destructor
  ~ComponentHeader() {
    delete [] info;
  }

  //! Write binary header to output stream
  bool write(ostream *out);

  //! Read binary header from input stream
  bool read(istream *in);

};

#endif
