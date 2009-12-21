// -*- C++ -*-

#ifndef _BarrierWrapper_H
#define _BarrierWrapper_H

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

#include <mpi.h>

class BarrierWrapper
{
private:
  char *buffer, *bufferT;
  bool check_label;
  MPI_Comm comm;
  int commsize, localid;
  static const int cbufsz = 128;

public:

  BarrierWrapper(MPI_Comm comm, bool label=false);
  BarrierWrapper(const BarrierWrapper& inst);
  ~BarrierWrapper();

  void operator()(const string &label);
};

#endif // _BarrierWrapper_H
