#ifndef _OutFrac_H
#define _OutFrac_H

#include <Output.H>
#include <Component.H>

/** Dump mass fractions at each interval
    
    @param filename of the fraction file

    @param name of the component to dump
  
    @param nint is the frequency between file updates 

    @param frac(n) are the quantiles
*/
class OutFrac : public Output
{

private:

  string filename;
  Component *tcomp;
  int numQuant;
  vector<double> Quant;

  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  OutFrac(const YAML::Node& conf);

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
