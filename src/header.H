#ifndef header_H
#define header_H

#include <iostream>
#include <iomanip>
#include <sstream>

#include <boost/shared_array.hpp>

#include "localmpi.h"

using namespace std;

//! Describe the phase space dump
class MasterHeader {
public:
  //! Current time for phase space
  double time;
  //! Number of particles in entire phase space
  int ntot;
  //! Number of individual components
  int ncomp;

  friend std::ostream& operator<<(std::ostream& os, const MasterHeader& p) 
  {
    os << std::setw(14) << std::right << "Time" << " : "
       << std::setw(14) << std::left  << p.time << std::endl
       << std::setw(14) << std::right << "Ntot" << " : "
       << std::setw(14) << std::left  << p.ntot << std::endl
       << std::setw(14) << std::right << "Ncomp" << " : "
       << std::setw(14) << std::left  << p.ncomp << std::endl;
    return os;
  }

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
  boost::shared_array<char> info;
  
  //! Default size of info string (if not specified explicitly or read)
  static int defaultInfoSize;

  //! Default constructor: uses defaultInfoSize by default
  ComponentHeader() 
  {
    ninfochar = defaultInfoSize;
    info = boost::shared_array<char>(new char [ninfochar]);
  }

  //! Construct component with info field of length n
  ComponentHeader(int n)
  {
    ninfochar = n;
    info = boost::shared_array<char>(new char [ninfochar]);
  }

  //! Write binary header to output stream
  bool write(ostream *out);

  //! Write binary header to MPI
  bool write_mpi(MPI_File& out, MPI_Offset& offset);

  //! Read binary header from input stream
  bool read(istream *in);

  //! Get header size
  int getSize()
  {
    return sizeof(int)*4 + ninfochar;
  }

  friend std::ostream& operator<<( std::ostream& os, const ComponentHeader& p ) 
  {
    os << std::setw(14) << std::right << "Nbod" << " : "
       << std::setw(14) << std::left  << p.nbod << std::endl
       << std::setw(14) << std::right << "N(int)" << " : "
       << std::setw(14) << std::left  << p.niatr << std::endl
       << std::setw(14) << std::right << "N(real)" << " : "
       << std::setw(14) << std::left  << p.ndatr << std::endl
       << std::setw(14) << std::right << "Info" << " : "
       << std::setw(60) << std::right << p.info.get() << std::endl;
    
    return os;
  }

};

#endif
