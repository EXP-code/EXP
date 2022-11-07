#ifndef _OutAscii_H
#define _OutAscii_H

#include <Component.H>
#include <OutAscii.H>

/** Write phase-space dumps for a single component 
    at regular intervals in ascii format.  Each line contains:
    mass, pos[3], vel[3], acc[3], pot, potext, ival[niattr], dval[ndattr]

    Each %dump is written into a new file label as <code>filename.n</code>
    where is begins at <code>nbeg</code> and incremented by 1 after each
    file is written.

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
    @param nintsub is the substep number to perform outputs
    @param nbeg is suffix of the first phase space %dump
    @param name of the desired component
    @param accel set to 1 means print out accelerations
*/
class OutAscii : public Output
{

private:

  string filename;
  string name;
  int nbeg;
  bool accel;

  Component *c0;

  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  OutAscii(const YAML::Node& conf);

  //! Provided by derived class to generate some output
  /*!
    \param nstep is the current time step used to decide whether or not
    to %dump
    \param last should be true on final step to force phase space %dump
    indepentently of whether or not the frequency criterion is met
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
