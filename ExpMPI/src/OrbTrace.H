#ifndef _OutTrace_H
#define _OutTrace_H

#include <OrbTrace.H>

/** Log norb orbits at each interval
    
    @param norb is the number of orbits per node to follow
  
    @param nskip is the interval between orbits beginning with nskip.
    If nskip is unspecified, <code>nskip=nbodies/norb</code>.
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
