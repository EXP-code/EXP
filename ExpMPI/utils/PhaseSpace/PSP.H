#ifndef _PSP_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <list>
#include <vector>
#include <string>
#include <cmath>

#include <boost/shared_ptr.hpp>

#include <header.H>
#include <StringTok.H>

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
  
template <typename real>
class PParticle 
{
public:
  real mass;
  real pos[3];
  real vel[3];
  real phi;
  vector<int>    iatr;
  vector<real>   datr;
  unsigned long  indx;

  void read(std::istream* in, int pcount, list<PSPstanza>::iterator spos) 
  {
    if (spos->index_size) 
      in->read((char *)&indx, sizeof(unsigned long));
    else
      indx = pcount;

    in->read((char *)&mass, sizeof(real));
    for (int i=0; i<3; i++) in->read((char *)&pos[i], sizeof(real));
    for (int i=0; i<3; i++) in->read((char *)&vel[i], sizeof(real));
    in->read((char *)&phi, sizeof(real));
    for (int i=0; i<spos->niatr; i++) 
      in->read((char *)&iatr[i], sizeof(int));
    for (int i=0; i<spos->ndatr; i++) 
      in->read((char *)&datr[i], sizeof(real));
  }

};

class SParticle
{
private:
  unsigned s;

public:
  boost::shared_ptr< PParticle<float>  > f;
  boost::shared_ptr< PParticle<double> > d;

  SParticle() : s(0) {}

  double mass() throw(std::string) {
    if (f.get()) return f->mass;
    if (d.get()) return d->mass;
    throw std::string("Particle not initialized");
    return 0.0;
  }

  double pos(int i) throw(std::string) {
    if (f.get()) return f->pos[i];
    if (d.get()) return d->pos[i];
    throw std::string("Particle not initialized");
    return 0.0;
  }

  double vel(int i) throw(std::string) {
    if (f.get()) return f->vel[i];
    if (d.get()) return d->vel[i];
    throw std::string("Particle not initialized");
    return 0.0;
  }

  double iatr(int i) throw(std::string) {
    if (f.get()) {
      if (i<f->iatr.size()) return f->iatr[i];
      else {
	throw std::string("Invalid index");
	return 0;
      }
    }
    if (d.get()) {
      if (i<d->iatr.size()) return d->iatr[i];
      else {
	throw std::string("Invalid index");
	return 0;
      }
    }
    throw std::string("Particle not initialized");
    return 0;
  }

  double datr(int i) throw(std::string) {
    if (f.get()) {
      if (i<f->datr.size()) return f->datr[i];
      else {
	throw std::string("Invalid index");
	return 0.0;
      }
    }
    if (d.get()) {
      if (i<d->datr.size()) return d->datr[i];
      else {
	throw std::string("Invalid index");
	return 0.0;
      }
    }
    throw std::string("Particle not initialized");
    return 0.0;
  }

  double phi() throw(std::string) {
    if (f.get()) return f->phi;
    if (d.get()) return d->phi;
    throw std::string("Particle not initialized");
    return 0.0;
  }

  int indx() throw(std::string) {
    if (f.get()) return f->indx;
    if (d.get()) return d->indx;
    throw std::string("Particle not initialized");
    return 0;
  }

  void setSize(unsigned rsize) 
  {
    if (s == rsize) return;
    if (s == sizeof(float) ) delete f.get();
    if (s == sizeof(double)) delete d.get();

    s = rsize;
    if (s == sizeof(float))
      f = boost::shared_ptr< PParticle<float>  >( new PParticle<float>() );
    else 
      d = boost::shared_ptr< PParticle<double> >( new PParticle<double>() );
  }

  void read(std::istream* in, unsigned rsize, int pcount, 
	    list<PSPstanza>::iterator spos) 
  {
    setSize(rsize);
    if (s == sizeof(float)) 
      f->read(in, pcount, spos);
    else
      d->read(in, pcount, spos);
  }

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
  int pcount, rsize;

  //! Temporaries for stanza statistics
  float mtot;
  std::vector<float> pmin, pmed, pmax;
  std::vector<float> vmin, vmed, vmax;
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

  SParticle *GetParticle (istream* in);
  SParticle *NextParticle(istream* in);
};

extern string trimLeft(const string);
extern string trimRight(const string);

#endif
