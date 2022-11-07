#ifndef _OutCalbr_H
#define _OutCalbr_H

#include <Output.H>
#include <Component.H>

/** Compute changes in energy and angular momentum per orbit 
    and provide rms changes as a function of energy
    
    @param filename of the output file

    @param name of the component to dump
  
    @param nint is the frequency between file updates 

    @param N are the number of energy bins

    The routine assumes that the input file is defined with four double
    attributes to hold the energy and angular momentum vector.
*/
class OutCalbr : public Output
{

private:

  string filename;
  Component *tcomp;
  int num;
  vector<double> Ec, deltaE, deltaLx, deltaLy, deltaLz;
  vector<double> deltaE1, deltaLx1, deltaLy1, deltaLz1;
  vector<unsigned> ncnt, ncnt1;
  double Emin, Emax, dE;

  void set_energies(void);
  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  OutCalbr(const YAML::Node& conf);

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
