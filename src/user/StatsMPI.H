
#ifndef StatsMPI_H
#define StatsMPI_H

#include "localmpi.h"
#include "Species.H"
#include <array>
#include <map>

class StatsMPI
{
private:
  std::map<speciesKey, double> cntT;
  std::map<speciesKey, double> curT;

  std::map<speciesKey, double> ccnt;
  std::map<speciesKey, double> vsum, vmin, vmax;

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
  void add(const speciesKey& k, double wgt, double val)
  { 
    if (cntT.find(k) == cntT.end()) {
      cntT[k]  = wgt;
      curT[k]  = val;
      vmin[k]  = val;
      vmax[k]  = val;
    } else {
      cntT[k] += wgt;
      curT[k] += val;
      vmin[k]  = std::min<double>(vmin[k], val);
      vmax[k]  = std::max<double>(vmax[k], val);
    }
  }

  typedef std::array<double, 4> SpStat;
  typedef std::map<speciesKey, SpStat> Return;

  //! Access values
  Return stats()
  {
    share();
    Return ret;
    for (auto v : ccnt) {
      speciesKey k = v.first;
      ret[k] = {ccnt[k], vsum[k], vmin[k], vmax[k]};
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
