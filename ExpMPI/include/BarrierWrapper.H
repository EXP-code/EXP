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

  //! Buffer size for checking working labels
  static int    cbufsz;

  //! Set to true while control is in the wrapper
  static bool   inOper;

  //! The label supplied by the caller
  static string lbOper;

  //! The source file containing the call
  static string flOper;

  //! The source file line number with the call
  static int    lnOper;

  //! Constructor
  BarrierWrapper(MPI_Comm comm, bool label=false);

  //! Copy constructor
  BarrierWrapper(const BarrierWrapper& inst);

  //! Destructor
  ~BarrierWrapper();

  void operator()(const string &label, const char* file, const int line);
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
