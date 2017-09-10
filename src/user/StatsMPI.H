
#ifndef StatsMPI_H
#define StatsMPI_H

#include "localmpi.h"
#include "Species.H"
#include <array>
#include <map>

class StatsMPI
{
private:
  std::map<speciesKey, unsigned long> cntT;
  std::map<speciesKey, double       > curT;

  std::map<speciesKey, unsigned long> ccnt;
  std::map<speciesKey, double       > vsum, vmin, vmax;

  void share()
  {
    for (auto s : cntT) {
      speciesKey k = s.first;
      if (ccnt.find(k) == ccnt.end()) {
	ccnt[k] = cntT[k];
	vsum[k] = curT[k];
      } else {
	ccnt[k] += cntT[k];
	vsum[k] += curT[k];
      }
    }

    cntT.clear();
    curT.clear();
  }

public:

  //! Constructor
  StatsMPI() {}
  //! Add a value
  void add(const speciesKey& k, double val)
  { 
    if (cntT.find(k) == cntT.end()) {
      cntT[k]  = 1;
      curT[k]  = val;
      vmin[k]  = val;
      vmax[k]  = val;
    } else {
      cntT[k]++;
      curT[k] += val;
      vmin[k]  = std::min<double>(vmin[k], val);
      vmax[k]  = std::max<double>(vmax[k], val);
    }
  }

  typedef std::array<double, 3> StatType;
  typedef std::pair<unsigned long, StatType> SpStat;
  typedef std::map<speciesKey, SpStat> Return;

  //! Access values
  Return stats()
  {
    share();
    Return ret;
    for (auto v : ccnt) {
      speciesKey k = v.first;
      ret[k] = {ccnt[k], {vsum[k], vmin[k], vmax[k]}};
    }
    return ret;
  }

  //! Reset values
  void clear()
  {
    ccnt.clear();
    vsum.clear();
    vmin.clear();
    vmax.clear();
  }

};

#endif
