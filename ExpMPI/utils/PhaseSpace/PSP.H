#ifndef _PSP_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <list>
#include <vector>
#include <string>
#include <cmath>

#include <header.H>

struct PSPstanza 
{
  streampos pos, pspos;
  string name;
  string id;
  string cparam;
  string fparam;
  string ttype;
  int nbod;
  int niatr;
  int ndatr;
  size_t index_size;
};
  
class SParticle 
{
public:
  double mass;
  double pos[3];
  double vel[3];
  double phi;
  vector<int>    iatr;
  vector<double> datr;
  unsigned long  indx;
};

class Dump 
{
public:
  streampos pos;
  MasterHeader header;
  list<PSPstanza> stanzas;
  
  int ngas, ndark, nstar, ntot;
  list<PSPstanza> gas, dark, star;
  
  Dump () : ngas(0), ndark(0), nstar(0), ntot(0) {}

  //! Add replace cparam field
  void add_cparam(void);
};

class PSPDump 
{

private:

  list<Dump> dumps;
  bool TIPSY, VERBOSE;
  
  list<Dump>::iterator sdump;
  Dump *fid;

  list<PSPstanza>::iterator spos;
  PSPstanza *cur;

  SParticle part;
  int pcount;

  //! Temporaries for stanza statistics
  float mtot;
  vector<float> pmin, pmed, pmax;
  vector<float> vmin, vmed, vmax;
  void ComputeStats(istream *in);

public:

  PSPDump(ifstream *in, bool tipsy=false, bool verbose=false);

  double SetTime(double time);
  double CurrentTime() {return fid->header.time;}
  void PrintSummaryCurrent(ifstream *in, ostream &out, bool stats=false, bool timeonly=false);
  void PrintSummary       (ifstream *in, ostream &out, bool stats=false, bool timeonly=false);

  Dump *CurrentDump() {return fid;}

  Dump *GetDump();
  Dump *NextDump();

  PSPstanza *GetStanza();
  PSPstanza *NextStanza();

  PSPstanza *GetGas();
  PSPstanza *NextGas();

  PSPstanza *GetDark();
  PSPstanza *NextDark();

  PSPstanza *GetStar();
  PSPstanza *NextStar();

  SParticle *GetParticle(istream* in);
  SParticle *NextParticle(istream* in);
};

#endif
