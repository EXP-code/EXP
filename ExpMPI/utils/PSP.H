#ifndef _PSP_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <list>
#include <vector>
#include <string>

#include <header.H>

struct Stanza {

  streampos pos, pspos;
  string name;
  string id;
  string param;
  string ttype;
  int nbod;
  int niatr;
  int ndatr;
};
  
class SParticle 
{
public:
  double mass;
  double pos[3];
  double vel[3];
  double phi;
  vector<int> iatr;
  vector<double> datr;
};

class Dump 
{
public:
  streampos pos;
  MasterHeader header;
  list<Stanza> stanzas;
  
  int ngas, ndark, nstar, ntot;
  list<Stanza> gas, dark, star;
  
  Dump () : ngas(0), ndark(0), nstar(0), ntot(0) {}
};

class PSPDump 
{

private:

  list<Dump> dumps;
  bool TIPSY;
  
  list<Dump>::iterator sdump;
  Dump *fid;

  list<Stanza>::iterator spos;
  Stanza *cur;

  SParticle part;
  int pcount;

public:

  PSPDump(ifstream *in, bool tipsy=false);

  double SetTime(double time);
  double CurrentTime() {return fid->header.time;}
  void PrintSummaryCurrent(ostream &out);
  void PrintSummary(ostream &out);

  Dump *CurrentDump() {return fid;}

  Dump *GetDump();
  Dump *NextDump();

  Stanza *GetStanza();
  Stanza *NextStanza();

  Stanza *GetGas();
  Stanza *NextGas();

  Stanza *GetDark();
  Stanza *NextDark();

  Stanza *GetStar();
  Stanza *NextStar();

  SParticle *GetParticle(istream* in);
  SParticle *NextParticle(istream* in);
};

#endif
