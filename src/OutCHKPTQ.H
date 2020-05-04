#ifndef _OUTCHKPTQ_H
#define _OUTCHKPTQ_H

#include <Output.H>

/** Writes a checkpoint file at regular intervals from each node in
    component pieces.  These pieces may be reassembled from the info
    in an "instruction" file.

    Previously written checkpoint file will be rename to
    \<filename\>.bak

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
    @param mpio set to true uses MPI-IO output with arbitrarily 
    sequenced particles
    @param nagg is the number of MPI-IO aggregators
*/
class OutCHKPTQ : public Output
{

private:

  std::string filename, nagg;
  bool timer, mpio;

  void initialize(void);

public:

  //! Constructor
  OutCHKPTQ(const YAML::Node& conf);

  //! Provided by derived class to generate some output
  /*!
    \param nstep is the current time step used to decide whether or not
    to dump
    \param last should be true on final step to force phase space dump
    indepentently of whether or not the frequency criterion is met
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
