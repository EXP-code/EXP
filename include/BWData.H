// -*- C++ -*-

#ifndef _BWData_H
#define _BWData_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <set>

using namespace std;

class BWData 
{
public:
  //! Number of seconds before first alarm
  static int dt1;

  //! Number of seconds before all sebsequent alarms
  static int dt2;

  typedef std::pair<time_t, int> Oelem;
  std::set<Oelem> owner;
  time_t first, expire;
  std::vector<bool> nd;
  unsigned count;

  BWData() {}

  BWData(int localid, int commsize) 
  {
    first       = time(0);
    expire      = first + dt1;
    count       = 1;
    nd          = std::vector<bool>(commsize);
    nd[localid] = true;

    owner.insert(Oelem(first, localid));
  }

  BWData(int localid, int commsize, time_t ctm) 
  {
    first       = ctm;
    expire      = first + dt1;
    count       = 1;
    nd          = std::vector<bool>(commsize);
    nd[localid] = true;

    owner.insert(Oelem(first, localid));
  }

  BWData(const BWData &p)
  {
    first  = p.first;
    expire = p.expire;
    nd     = p.nd;
    count  = p.count;
    owner  = p.owner;
  }

  void Add(int localid, time_t ctm) { owner.insert(Oelem(ctm, localid)); }

  int CTime() { return owner.begin()->first;  }
  int Owner() { return owner.begin()->second; }
};

#endif
