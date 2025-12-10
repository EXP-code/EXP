#ifndef _OutRelaxation_H
#define _OutRelaxation_H

#include "Output.H"

//! Output routine for generateRelaxation ("external force") module
class OutRelaxation : public Output
{

private:

  std::string fname;
  std::string suffix;

  int epos;
  void initialize(void);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  OutRelaxation(const YAML::Node& conf);
  //! Destructor
  virtual ~OutRelaxation() {};

  //! Do some output
  /*!
    \param nstep is the current time step index
    \param mstep is the current multistep level to decide whether or not to dump multisteps
    \param final should be true on last step to force output
  */
  void Run(int nstep, int mstep, bool final=false);

};

#endif
