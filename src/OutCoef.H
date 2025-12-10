#ifndef _OutCoef_H
#define _OutCoef_H

#include "OutCoef.H"

/** Dump coefficients at each interval
    
    @param filename of the coefficient file

    @param name of the component to dump
  
    @param nint is the frequency between file updates 

    @param native set to true uses old-style native coefficient format
*/
class OutCoef : public Output
{

private:

  std::string filename;
  double prev = -std::numeric_limits<double>::max();
  Component *tcomp;
  bool native;

  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  OutCoef(const YAML::Node& conf);

  //! Generate the output
  /*!
    \param nstep is the current time step used to decide whether or not
    to dump
    \param last should be true on final step to force phase space dump
    indepentently of whether or not the frequency criterion is met
  */
  void Run(int nstep, int mstep, bool last);

};

#endif
