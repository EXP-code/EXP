// -*- C++ -*-

#ifndef _BWHelper_H
#define _BWHelper_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <limits>
#include <set>

using namespace std;

class bw_struct 
{
public:
  //! Number of seconds before first alarm
  static int dt1;

  //! Number of seconds before all sebsequent alarms
  static int dt2;

  time_t first, expire;
  std::vector<bool> nd;
  unsigned count;
  int owner;

  bw_struct(int localid, int commsize) 
  {
    first       = time(0);
    expire      = first + dt1;
    owner       = localid;
    count       = 1;
    nd          = std::vector<bool>(commsize);
    nd[localid] = true;
  }

  bw_struct(int owner, time_t first)
  {
    this->owner  = owner;
    this->first  = first;
    this->expire = INT_MAX;
  }

};

#endif
