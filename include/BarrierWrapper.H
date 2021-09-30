// -*- C++ -*-

#ifndef _BarrierWrapper_H
#define _BarrierWrapper_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <set>

using namespace std;

#include <memory>

#include <Timer.H>
#include <mpi.h>


typedef std::shared_ptr<char> CharPtr;
typedef std::shared_ptr<MPI_Request> ReqPtr;

class Info 
{
private:

  //! Request handle
  std::list<ReqPtr> req;

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

  //! Construct from MPI-passed char buffer
  Info(CharPtr p);

  //! Full constructor
  Info(int own, time_t ctm, const std::string& s);

  //! Send the Info to all processes in the communicator
  void sendInfo(int myid, MPI_Comm comm, int commsize);

};

typedef std::shared_ptr<Info> InfoPtr;


class BWData 
{
public:
  //! Number of seconds before first alarm
  static int dt1;

  //! Number of seconds before all subsequent alarms
  static int dt2;

  //@{
  //! Times: initial and next expiration
  time_t first, expire;
  //@}

  //! Array of received entries
  std::vector<bool> nd;

  //! Number of alarm expirations
  unsigned nexpiry;

  //! Number of entries registered
  unsigned count;

  //! Cache Info shared pointers; storage reclaimed on destruction
  std::list<InfoPtr> info;

  //! Full constructor
  BWData(InfoPtr& p, int commsize);

  //! Register node in existing entry
  void Add(InfoPtr& p);

  //! Return the valid node with the lowest rank
  int Owner();
};


typedef std::shared_ptr<BWData> BWPtr;

class BarrierWrapper
{
  friend class Info;

private:

  typedef std::pair<CharPtr, ReqPtr> RecvElem;
  std::list<RecvElem> req;
  unsigned queued;
  unsigned nrecv;

  char *buffer, *bufferT;
  int commsize, localid;
  bool check_label;
  MPI_Comm comm;
  Timer timer;
  bool onoff;

  void light_operator(const string &label, const char* file, const int line);
  void heavy_operator(const string &label, const char* file, const int line);
  static bool light;

  std::map<std::string, BWPtr> pending;

  void updateMap(InfoPtr p);
  void syncTest(const std::string& mesg, const std::string& label);
  void listReport(const char*, std::map<std::string, BWPtr>::iterator);
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
    double ret = timer.getTime();
    timer.reset();
    return ret;
  };

  std::vector<int> getMissing(BWPtr p);

};



#endif // _BarrierWrapper_H
