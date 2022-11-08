#ifndef _BADNELL_H
#define _BADNELL_H

#include <utility>
#include <string>
#include <vector>
#include <map>

#include <memory>

class atomicData;
class Ion;

/** 
   Reads Nigel Badnell's cross-section data files

   Enter the (Z, C) pairs in ionQ database (see constructor
   implementation) and put the files in the RR and DR directories.
   The data directory is ./ by default and can be set by the
   environment variable ATOMIC_DATA_PATH.

   For ions not in the database, it will use the hydrogenic RR cross
   section.  This method can be used with RECOMB_RATIO in CollideIon.
   In this case, the ADAS cross sections will be used where they exist
   and the hydrogenic cross section with the fudge factor will be used
   for others.
*/
class BadnellData
{
public:

  //! Pair of unsigned shorts for indexing (Z, C) index pairs
  using lQ = std::pair<unsigned short, unsigned short>;

  //! Rydberg energy in eV
  static constexpr double RydtoeV = 13.60569253;

private:

  atomicData *ad;

  //! Ion mnemonics for file names
  std::map<int, std::string> ions;

  //! Ions with file data in the RR and DR directories in the datapath
  std::vector<lQ> ionQ;

  //! Parse directory
  std::tuple<std::map<int, std::string>, std::vector<lQ>> walkDirectory();

public:

  //! Location of Nigel Badnell's adf files in RR and DR
  //! subdirectories for the radiative and dielectronic recombination
  //! cross sections.  Environment variable ATOMIC_DATA_PATH can be
  //! used to set this variable.
  static std::string datapath;

  //! Reweight cross section by mean energy (default: false).  Doesn't
  //! matter much for small bins.
  static bool reweight;
  
  //! Coupling choice (default intermediate coupling, "ic")
  static std::string coupling;

  //! Store the cross-section data table for each ion
  struct BadnellRec 
  {
    std::vector<double> E_rr, X_rr, E_dr, X_dr;

    void sync(int myid=0);
  };
  
  using bdPtr = std::shared_ptr<BadnellRec>;

  std::map<lQ, bdPtr> data;

  //! Constructor
  BadnellData();

  //! Initialize the data from the data file
  void initialize(atomicData* ad);

  //! Calculates the summed RR and DR recombination cross section
  double cross(const lQ& Q, double EeV);

};

#endif
