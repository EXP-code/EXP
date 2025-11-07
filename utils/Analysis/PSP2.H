#ifndef _PSP2_H

#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <memory>
#include <string>
#include <cmath>
#include <list>

#include <memory>

#include "StringTok.H"
#include "header.H"

class PSPstanza 
{
public:
  ComponentHeader comp;
  std::string name;
  std::string id;
  std::string cparam;
  std::string fparam;
  size_t index_size, r_size;

  streampos pos, pspos;
  std::vector<std::string> nparts;

  bool operator==(const PSPstanza& x) const
  {
    if(name.compare(x.name) == 0)
      return true;
    else
      return false;
  }

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

  void read(std::istream& in, int pcount, list<PSPstanza>::iterator spos) 
  {
    // Sequence value
    // --------------
    if (spos->index_size) 
      in.read((char *)&indx, sizeof(unsigned long));
    else
      indx = pcount;

    // Initialize fields
    // -----------------
    if (spos->comp.niatr) {
      if (iatr.size() != spos->comp.niatr)
	iatr = vector<int>(spos->comp.niatr);
    }
    if (spos->comp.ndatr) {
      if (datr.size() != spos->comp.ndatr)
	datr = vector<real>(spos->comp.ndatr);
    }


    // Read the rest
    // -------------
    in.read((char *)&mass, sizeof(real));
    for (int i=0; i<3; i++) in.read((char *)&pos[i], sizeof(real));
    for (int i=0; i<3; i++) in.read((char *)&vel[i], sizeof(real));
    in.read((char *)&phi, sizeof(real));
    for (int i=0; i<spos->comp.niatr; i++) 
      in.read((char *)&iatr[i], sizeof(int));
    for (int i=0; i<spos->comp.ndatr; i++) 
      in.read((char *)&datr[i], sizeof(real));
  }

};


class SParticle
{
private:
  unsigned s;

public:
  std::shared_ptr< PParticle<float>  > f;
  std::shared_ptr< PParticle<double> > d;

  SParticle() : s(0) {}

  double mass()
  {
    if (f.get()) return f->mass;
    if (d.get()) return d->mass;
    throw std::runtime_error("Particle not initialized");
    return 0.0;
  }

  double pos(int i)
  {
    if (f.get()) return f->pos[i];
    if (d.get()) return d->pos[i];
    throw std::runtime_error("Particle not initialized");
    return 0.0;
  }

  double vel(int i)
  {
    if (f.get()) return f->vel[i];
    if (d.get()) return d->vel[i];
    throw std::runtime_error("Particle not initialized");
    return 0.0;
  }

  int niatr()
  {
    if (f.get()) return f->iatr.size();
    if (d.get()) return d->iatr.size();
    throw std::runtime_error("Particle not initialized");
    return 0;
  }

  int ndatr()
  {
    if (f.get()) return f->datr.size();
    if (d.get()) return d->datr.size();
    throw std::runtime_error("Particle not initialized");
    return 0;
  }

  int iatr(int i)
  {
    if (f.get()) {
      if (i<f->iatr.size()) return f->iatr[i];
      else {
	throw std::runtime_error("Invalid index");
	return 0;
      }
    }
    if (d.get()) {
      if (i<d->iatr.size()) return d->iatr[i];
      else {
	throw std::runtime_error("Invalid index");
	return 0;
      }
    }
    throw std::runtime_error("Particle not initialized");
    return 0;
  }

  double datr(int i)
  {
    if (f.get()) {
      if (i<f->datr.size()) return f->datr[i];
      else {
	throw std::runtime_error("Invalid index");
	return 0.0;
      }
    }
    if (d.get()) {
      if (i<d->datr.size()) return d->datr[i];
      else {
	throw std::runtime_error("Invalid index");
	return 0.0;
      }
    }
    throw std::runtime_error("Particle not initialized");
    return 0.0;
  }

  double phi()
  {
    if (f.get()) return f->phi;
    if (d.get()) return d->phi;
    throw std::runtime_error("Particle not initialized");
    return 0.0;
  }

  int indx()
  {
    if (f.get()) return f->indx;
    if (d.get()) return d->indx;
    throw std::runtime_error("Particle not initialized");
    return 0;
  }

  void setSize(unsigned rsize) 
  {
    if (s == rsize) return;
    if (s == sizeof(float) ) delete f.get();
    if (s == sizeof(double)) delete d.get();

    s = rsize;
    if (s == sizeof(float))
      f = std::shared_ptr< PParticle<float>  >( new PParticle<float>() );
    else 
      d = std::shared_ptr< PParticle<double> >( new PParticle<double>() );
  }

  void read(std::istream& in, unsigned rsize, int pcount, 
	    list<PSPstanza>::iterator spos) 
  {
    setSize(rsize);
    if (s == sizeof(float)) 
      f->read(in, pcount, spos);
    else
      d->read(in, pcount, spos);
  }

  void write(std::ostream& out, bool real4, size_t isiz);
};

class PSP
{

protected:

  const static unsigned long magic = 0xadbfabc0;
  const static unsigned long mmask = 0xf;
  const static unsigned long nmask = ~mmask;

  MasterHeader header;

  bool VERBOSE;
  
  list<PSPstanza> stanzas;
  list<PSPstanza>::iterator spos;
  PSPstanza *cur;

  SParticle part;
  int pcount;

  std::ifstream in;

  //! Temporaries for stanza statistics
  float mtot;
  std::vector<float> pmin, pmed, pmax;
  std::vector<float> vmin, vmed, vmax;
  void ComputeStats();
  void write_binary(std::ostream& out,
		    list<PSPstanza>::iterator its, bool real4);

  //! Append '/' if necessary
  void check_dirname();
  std::string new_dir;

  //! Stream init
  void init()
  {
    // Prepare <in> to throw if failbit gets set
    in.exceptions(std::ifstream::failbit | std::ifstream::badbit);
  }

public:

  //! Default constructor
  PSP(bool verbose) : VERBOSE(verbose), new_dir("./") { init(); }

  //! Constructor specifying current working data directory
  PSP(bool verbose, std::string new_dir) : VERBOSE(verbose), new_dir(new_dir)
  { init(); check_dirname(); }

  //! Destructor
  virtual ~PSP() { if (in.is_open()) in.close(); }

  //! Get snapshot time
  double CurrentTime()    {return header.time;}
  
  //! Print summary phase-space info
  void   PrintSummary     (ostream &out, bool stats=false, bool timeonly=false);

  //! Set stanza to name
  virtual PSPstanza *GetNamed(const std::string& name)
  {
    PSPstanza x;
    x.name = name;
    auto it = std::find(stanzas.begin(), stanzas.end(), x);
    if (it == stanzas.end()) return 0;
    spos = it;
    cur  = &(*spos);
    return cur;
  }

  //! Reset to beginning of stanza list and return null if no stanzas.
  virtual PSPstanza *GetStanza();

  //! Get next stanza in list and return null at end;
  virtual PSPstanza *NextStanza();

  //! Reset to beginning of particles for this component
  virtual SParticle *GetParticle () = 0;

  //! Get next particle and return null at end.
  virtual SParticle *NextParticle() = 0;

  //! Write a new PSP file
  void writePSP(std::ostream& out,  bool real4);

};



/**
   Class to access a full PSP file (OUT)
 */
class PSPout : public PSP
{
private:
  streampos pos;

public:
  PSPout(const std::string& file, bool verbose=false);

  //! Destructor
  virtual ~PSPout() {}

  virtual SParticle* GetParticle ();
  virtual SParticle* NextParticle();
};

/**
   Class to access a SPLIT PSP file (SPL)
 */
class PSPspl : public PSP
{
private:
  unsigned int N;
  unsigned int fcount;
  std::vector<std::string>::iterator fit;

  //! Open next file part
  void openNextBlob();

public:
  //! Constuctors
  PSPspl(const std::string& file, const std::string dir = "", bool verbose=false);

  //! Destructor
  virtual ~PSPspl() {}

  virtual SParticle* GetParticle ();
  virtual SParticle* NextParticle();

};

extern std::string trimLeft(const string);
extern std::string trimRight(const string);

typedef std::shared_ptr<PSP> PSPptr;

#endif
