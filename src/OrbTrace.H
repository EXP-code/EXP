#ifndef _OrbTrace_H
#define _OrbTrace_H

#include "OrbTrace.H"

/** Log norb orbits at each interval
    
    For best efficiency, particles to be traced should be 
    on the root process.

    @param norb is the number of orbits per node to follow
  
    @param nbeg is the offset in the bodies list for the first
    particle to trace.
    
    @param nskip is the interval between orbits beginning with nskip.
    If nskip is unspecified, <code>nskip=nbodies/norb</code>. 

    @param nint is the frequency between file updates 

    @param use_acc to output the acceleration

    @param use_pot to output the potential

    @param use_lev to output the multistep level

    @param orbitlist is the list of particle numbers to trace

    @param name of the component to trace
*/
class OrbTrace : public Output
{

private:

  std::string filename;
  std::string orbitlist;
  int nbeg;
  int nskip;
  int norb;
  bool use_acc;
  bool use_pot;
  bool use_lev;
  bool local;
  Component *tcomp;
  std::vector<int> orblist;
  std::vector<double> pbuf;
  int nbuf;
  int flags;

  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

  //! Last output time
  double prev = -std::numeric_limits<double>::max();

public:

  //! Constructor
  OrbTrace(const YAML::Node& conf);

  //! Generate the output
  /*!
    \param nstep is the current time step used to decide whether or not
    to dump
    \param mstep is the current multistep level to decide whether or not to dump multisteps
    \param last should be true on final step to force phase space dump
    indepentently of whether or not the frequency criterion is met
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
