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
  //! Pair of unsigned shorts for indexing (Z, C) index pairs
  using lQ = std::pair<unsigned short, unsigned short>;

private:

  chdata *ch;

  // Ion mnemonics for file names
  std::map<int, std::string> ions;

  // Members
  std::vector<lQ> ionQ;

public:

  //! Location of RR data files
  static std::string RRdir;

  //! Location of DR data files
  static std::string DRdir;

  //! Reweight cross section by mean energy
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
