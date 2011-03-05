// -*- C++ -*-

#ifndef _BarrierWrapper_H
#define _BarrierWrapper_H

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#include <Timer.h>
#include <mpi.h>

class BarrierWrapper
{
private:
  char *buffer, *bufferT;
  bool check_label;
  MPI_Comm comm;
  int commsize, localid;
  Timer timer;
  bool onoff;

public:

  static int cbufsz;

  BarrierWrapper(MPI_Comm comm, bool label=false);
  BarrierWrapper(const BarrierWrapper& inst);
  ~BarrierWrapper();

  void operator()(const string &label);
  void on()  { onoff = true;  }
  void off() { onoff = false; }

  double getTime() {
    // double ret = timer.getTime().getRealTime()*1.0e-6;
    double ret = timer.getTime()()*1.0e-6;
    timer.reset();
    return ret;
  };

};

#endif // _BarrierWrapper_H
