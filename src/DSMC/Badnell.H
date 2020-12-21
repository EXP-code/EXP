#ifndef _BADNELL_H
#define _BADNELL_H

#include <utility>
#include <string>
#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

class chdata;
class Ion;

/** 
   Reads Nigel Badnell's cross-section data files
*/
class BadnellData
{
public:

  //! Pair of unsigned shorts for indexing (Z, C) index pairs
  using lQ = std::pair<unsigned short, unsigned short>;

  //! Rydberg energy in eV
  static constexpr double RydtoeV = 13.60569253;

private:

  chdata *ch;

  // Ion mnemonics for file names
  std::map<int, std::string> ions;

  // Members
  std::vector<lQ> ionQ;

public:

  //! Location of Nigel Badnell's adf files in RR and DR
  //! subdirectories for the radiative and dielectronic recombination
  //! cross sections.  Environment variable ADF_DATA_PATH can be used
  //! to set this variable.
  static std::string datapath;

  //! Reweight cross section by mean energy (default: false).  Doesn't
  //! matter much for small bins.
  static bool reweight;
  
  struct BadnellRec 
  {
    std::vector<double> E_rr, X_rr, E_dr, X_dr;

    void sync(int myid=0);
  };
  
  using bdPtr = boost::shared_ptr<BadnellRec>;

  std::map<lQ, bdPtr> data;

  //! Constructor
  BadnellData();

  //! Initialize the data from the data file
  void initialize(chdata* ch);

  //! Calculates the summed recombination cross section
  double cross(const lQ& Q, double EeV);

};

#endif
