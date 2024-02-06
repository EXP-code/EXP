#ifndef _OutDiag_H
#define _OutDiag_H

#include <OutDiag.H>

//! Print field quantities for spherical and cylindrical bases
class OutDiag : public Output
{

private:

  std::string filename;
  std::list<Component*> lcomp;
  std::vector<std::string> names;
  double prev = -std::numeric_limits<double>::max();

  int NUM;
  double THETA, PHI, RMIN, RMAX;

  void initialize(void);
  void header(ostream& out);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  OutDiag(const YAML::Node& conf);

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
