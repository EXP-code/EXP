#ifndef _OutTrace_H
#define _OutTrace_H

#include <OrbTrace.H>

/** Log norb orbits at each interval
    
    @param norb is the number of orbits per node to follow
  
    @param nskip is the interval between orbits beginning with nskip.
    If nskip is unspecified, <code>nskip=nbodies/norb</code>.  The
    first particle traced will be offset in the bodies list by
    <code>nskip</code>.  If <code>nskip*norb</code> is larger than
    the number of bodies, only <code>nbodies/nskip</code> orbits
    will be traced.

    @param nfreq is the frequency between file updates 
*/
class OrbTrace : public Output
{

private:

  string filename;
  int nskip;
  int norb;

  void initialize(void);

public:

  //! Constructor
  OrbTrace(string& line);

  //! Generate the output
  /*!
    <param> nstep is the current time step used to decide whether or not
    to dump
    <param> last should be true on final step to force phase space dump
    indepentently of whether or not the frequency criterion is met
  */
  void Run(int nstep, bool last);

};

#endif
