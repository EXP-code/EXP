// -*- C++ -*-

#ifndef _BarrierWrapper_H
#define _BarrierWrapper_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

#include <boost/shared_ptr.hpp>

#include <BWData.H>
#include <Timer.h>
#include <mpi.h>


class Info;

class BarrierWrapper
{
  friend class Info;

private:

  char *buffer, *bufferT;
  int commsize, localid;
  bool check_label;
  MPI_Comm comm;
  Timer timer;
  bool onoff;

  void light_operator(const string &label, const char* file, const int line);
  void heavy_operator(const string &label, const char* file, const int line);
  static bool light;

  std::map<std::string, BWData> pending;

  void listReport(const char*, std::map<std::string, BWData>::iterator);
  void finalReport(std::string&);

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

  //@{
  //! Debugging output flag
  static bool   verbose;
  static bool   extra_verbose;
  static bool   debugging;
  //@}

  //! Loop wait in microseconds for heavy operator mode (default: 100)
  static int    loop_delay;

  //! Constructor
  BarrierWrapper(MPI_Comm comm, bool label=false);

  //! Copy constructor
  BarrierWrapper(const BarrierWrapper& inst);

  //! Destructor
  ~BarrierWrapper();

  void operator()(const string &label, const char* file, const int line)
  {
    if (light) light_operator(label, file, line);
    else       heavy_operator(label, file, line);
  }

  void on()  { onoff = true;  }
  void off() { onoff = false; }

  void setLightWeight() { light = true;  }
  void setHeavyWeight() { light = false; }

  double getTime() {
    // double ret = timer.getTime().getRealTime()*1.0e-6;
    double ret = timer.getTime()()*1.0e-6;
    timer.reset();
    return ret;
  };

};

typedef boost::shared_ptr<char> CharPtr;

class Info 
{
private:

  //! MPI buffer
  CharPtr blob;

  //! Helper to translate buffer
  void pack();

  //! Space for character string
  size_t charsz;

public:

  //! Buffer size (1 block)
  const static size_t bufsz = 512;

  //! Tag
  const static int tag = 232395;

  //! Owning process
  int own;

  //! Creation time
  time_t ctm;

  //! Label size
  unsigned siz;

  //! Label
  std::string s;

  //! Character buffer for string
  CharPtr c;

  //! Null constructor for assignment
  Info() {}

  Info(CharPtr p) {
    unpackInfo(p);
  }

  //! Full constructor
  Info(int own, time_t ctm, const std::string& s)
  {
    this->own = own;
    this->ctm = ctm;
    this->s   = s;
    this->siz = s.size() + 1;
    this->c   = CharPtr(new char [siz]);
    strncpy(c.get(), s.c_str(), siz);
  }

  //! Send the Info to all processes in the communicator
  void sendInfo(int myid, MPI_Comm comm, int commsize);

  //! Get the Info from the character blob
  void unpackInfo(CharPtr p);

};

typedef boost::shared_ptr<Info> InfoPtr;


#endif // _BarrierWrapper_H
