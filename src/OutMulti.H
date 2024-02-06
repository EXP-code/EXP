#ifndef _OutMulti_H
#define _OutMulti_H

#include <OutMulti.H>

/** Dump histogram of multistep levels to file in GNUPLOT format

    @param filename is the name of the output file
    @param nint is the number of steps between dumps
*/
class OutMulti : public Output
{

private:

  std::string filename;
  void initialize(void);

  double prev = -std::numeric_limits<double>::max();

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  OutMulti(const YAML::Node& conf);

  //! Provided by derived class to generate some output
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
