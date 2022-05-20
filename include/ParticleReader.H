#ifndef _ParticleReader_H
#define _ParticleReader_H

#include <unordered_map>
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

#include <mpi.h>

#include <config.h>

#include <StringTok.H>
#include <header.H>
#include <Particle.H>
#include <gadget.H>

#ifdef HAVE_XDR
#include <tipsy.H>
#endif

namespace PR
{

  //! Base class for reading particle phase space from any simulation
  class ParticleReader
  {
  protected:
    
    static std::vector<std::string> readerTypes;
    int numprocs, myid;
    
  public:
    
    //! Constructor: check for and set up MPI
    ParticleReader()
    {
      // Check for MPI
      int flag;
      MPI_Initialized(&flag);
      if (flag) {
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
      } else {
	numprocs = 1;
	myid     = 0;
      }
    }
    
    //! Select a particular particle type and reset the iterator
    virtual void SelectType(const std::string& type) = 0;
    
    //! Number of particles in the chosen type
    virtual unsigned long CurrentNumber() = 0;
    
    //! Return list of particle types
    virtual std::vector<std::string> GetTypes() = 0;
    
    //! Get current time
    virtual double CurrentTime() = 0;
    
    //! Reset to beginning of particles for this component
    virtual const Particle* firstParticle() = 0;
    
    //! Get the next particle
    virtual const Particle* nextParticle() = 0;
    
    //! Order file list into batches at single times
    static std::vector<std::vector<std::string>>
    parseFileList(const std::string& file, const std::string& delimit);

    //! Create a shared pointer to derived instance of the desired type and return it
    static std::shared_ptr<ParticleReader>
    createReader(const std::string& type, const std::vector<std::string>& file, int myid,
		 bool verbose=false);
  };
  
  class GadgetNative : public ParticleReader
  {
    gadget_header header;
    unsigned long totalCount;
    
    double time;
    std::vector<Particle> particles;
    bool _verbose;
    std::vector<std::string> _files;

    double mass[6];
    int npart[6], nptot[6];

    //! Current file
    std::vector<std::string>::iterator curfile;

    static std::vector<std::string> Ptypes;
    static std::unordered_map<std::string, int> findP;
    std::vector<std::string> Pfound;
    int ptype;
    
    unsigned pcount;
    void read_and_load();
    
    void getNumbers();
    void packParticle();
    bool nextFile();

  public:
    
    //! Constructor
    GadgetNative(const std::vector<std::string>& file, bool verbose=false);
    
    //! Select a particular particle type and reset the iterator
    void SelectType(const std::string& type) {
      auto it = findP.find(type);
      if (it != findP.end()) ptype = it->second;
      else {
	std::cerr << "GadgetNative error: no particle type <" << type
		  << ">" << std::endl;
	std::cerr << "Valid GadgetNative types are:";
	for (auto s : Ptypes) std::cerr << " " << s;
	std::cerr << std::endl;
      }

      curfile = _files.begin();	// Set to first file and open
      nextFile();
    }
    
    //! Number of particles in the chosen type
    unsigned long CurrentNumber() { return totalCount; }
    
    //! Return list of particle types
    std::vector<std::string> GetTypes() { return Pfound; }
    
    //! Get current time
    double CurrentTime() { return time; }
    
    //! Reset to beginning of particles for this component
    virtual const Particle* firstParticle();
    
    //! Get the next particle
    virtual const Particle* nextParticle();
    
  };
  
  
  class GadgetHDF5 : public ParticleReader
  {
  protected:
    unsigned long totalCount;
    
    double time;
    std::vector<Particle> particles;
    bool _verbose;
    std::vector<std::string> _files;

    double mass[6];
    int npart[6], nptot[6];
    
    
    //! Current file
    std::vector<std::string>::iterator curfile;

    static std::vector<std::string> Ptypes;
    static std::unordered_map<std::string, int> findP;
    std::vector<std::string> Pfound;
    int ptype;
    
    unsigned pcount;
    void read_and_load();
    
    void getNumbers();
    void packParticle();
    bool nextFile();

  public:
    
    //! Constructor
    GadgetHDF5(const std::vector<std::string>& file, bool verbose=false);
    
    //! Select a particular particle type and reset the iterator
    virtual void SelectType(const std::string& type) {
      auto it = findP.find(type);
      if (it != findP.end()) ptype = it->second;
      else {
	std::cerr << "GadgetHDF5 error: no particle type <" << type << ">" << std::endl;
	std::cerr << "Valid GadgetHDF5 types are:";
	for (auto s : Ptypes) std::cerr << " " << s;
	std::cerr << std::endl;
      }

      curfile = _files.begin();	// Set to first file and open
      nextFile();
    }
    
    //! Number of particles in the chosen type
    virtual unsigned long CurrentNumber() { return totalCount; }
    
    //! Return list of particle types
    virtual std::vector<std::string> GetTypes() { return Pfound; }
    
    //! Get current time
    virtual double CurrentTime() { return time; }
    
    //! Reset to beginning of particles for this component
    virtual const Particle* firstParticle();
    
    //! Get the next particle
    virtual const Particle* nextParticle();
    
  };
  
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
  class PParticle : public Particle
  {
  private:
    
    real _mass;
    real _pos[3];
    real _vel[3];
    real _phi;
    
    vector<real>   _datr;
    
  public:
    
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
	iattrib.resize(spos->comp.niatr);
      }
      if (spos->comp.ndatr) {
	_datr.resize(spos->comp.ndatr);
	dattrib.resize(spos->comp.ndatr);
      }
      
      
      // Read the rest
      // -------------
      in.read((char *)&_mass, sizeof(real));
      for (int i=0; i<3; i++) in.read((char *)&_pos[i], sizeof(real));
      for (int i=0; i<3; i++) in.read((char *)&_vel[i], sizeof(real));
      in.read((char *)&_phi, sizeof(real));
      for (int i=0; i<spos->comp.niatr; i++) 
	in.read((char *)&iattrib[i], sizeof(int));
      for (int i=0; i<spos->comp.ndatr; i++) 
	in.read((char *)&_datr[i], sizeof(real));
      
      // Copy real to double
      // -------------------
      mass = _mass;
      for (int i=0; i<3; i++) {
	pos[i] = _pos[i];
	vel[i] = _vel[i];
      }      
      
      dattrib.clear();
      for (auto v : _datr) dattrib.push_back(v);
      _datr.clear();
    }
    
    void skip(std::istream& in, int pcount, list<PSPstanza>::iterator spos) 
    {
      unsigned skipsize = 8*sizeof(real) +
	spos->comp.niatr*sizeof(int) +
	spos->comp.ndatr*sizeof(real);
      
      if (spos->index_size) skipsize += sizeof(unsigned long);
      
      in.seekg(skipsize, std::ios::cur);
    }
    
  };
  
  
  class PSP : public ParticleReader
  {
    
  protected:
    
    const static unsigned long magic = 0xadbfabc0;
    const static unsigned long mmask = 0xf;
    const static unsigned long nmask = ~mmask;
    
    MasterHeader header;
    
    bool VERBOSE;
    
    std::list<PSPstanza> stanzas;
    std::list<PSPstanza>::iterator spos;
    PSPstanza *cur;
    
    PParticle<float>   fpart;
    PParticle<double>  dpart;
    
    int pcount;
    
    std::ifstream in;
    
    //! Temporaries for stanza statistics
    float mtot;
    std::vector<double> pmin, pmed, pmax;
    std::vector<double> vmin, vmed, vmax;
    void ComputeStats();
    void write_binary(std::ostream& out,
		      list<PSPstanza>::iterator its, bool real4);
    
    //! Stream init
    void init()
    {
      // Prepare <in> to throw if failbit gets set
      in.exceptions(std::ifstream::failbit | std::ifstream::badbit);
      
    }
    
  public:
    
    //! Default constructor
    PSP(bool verbose) : VERBOSE(verbose) { init(); }
    
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
    
    std::vector<std::string> GetTypes()
    {
      std::vector<std::string> ret;
      for (auto v : stanzas) ret.push_back(v.name);
      return ret;
    }
    
    //! Set to type
    virtual void SelectType(const std::string& name)
    {
      if (GetNamed(name)) return;
      std::cout << "PSP error: no particle type <" << name << ">" << std::endl;
    }
    
    //! Number of particles in current type
    unsigned long CurrentNumber() { return cur->comp.nbod; }
    
    //! Reset to beginning of stanza list and return null if no stanzas.
    virtual PSPstanza *GetStanza();
    
    //! Get next stanza in list and return null at end;
    virtual PSPstanza *NextStanza();
    
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
    
    //! Constructor
    PSPout(const std::vector<std::string>& file, bool verbose=false);
    
    //! Destructor
    virtual ~PSPout() {}
    
    //@{
    //! Particle access
    virtual const Particle* firstParticle();
    virtual const Particle* nextParticle ();
    //@}
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
    PSPspl(const std::vector<std::string>& file, bool verbose=false);
    
    //! Destructor
    virtual ~PSPspl() {}
    
    //@{
    //! Particle access
    virtual const Particle* firstParticle ();
    virtual const Particle* nextParticle();
    //@}
    
  };
  
#ifdef HAVE_XDR
  /**
     Class to access a Tipsy file
  */
  class Tipsy : public ParticleReader
  {
  private:
    //! List of files
    std::vector<std::string> files;

    //! Current file
    std::vector<std::string>::iterator curfile;

    //! Pointer to Tipsy particle structures
    std::shared_ptr<TipsyXDR::TipsyFile> ps;

    //! Temporary for packing array
    Particle P;
    
    //! Current types
    std::vector<std::string> curTypes;
    
    //! Current name
    std::string curName;
    
    //! Totals
    unsigned long Ngas, Ndark, Nstar;

  protected:

    static std::vector<std::string> Ptypes;
    static std::unordered_map<std::string, int> findP;
    int ptype;

    unsigned pcount;

    void getNumbers();
    void packParticle();
    bool nextFile();

  public:
    
    //! Single-file constructor
    Tipsy(const std::string& file, bool verbose=false);
    
    //! Constructor
    Tipsy(const std::vector<std::string>& file, bool verbose=false);
    
    //! Destructor
    virtual ~Tipsy() {}
    
    //! Return list of particle types
    virtual std::vector<std::string> GetTypes() { return curTypes; }

    //! Select a particular particle type and reset the iterator
    virtual void SelectType(const std::string& type);

    //! Number of particles in the chosen type
    virtual unsigned long CurrentNumber();
    
    //! Get current time
    virtual double CurrentTime() { return ps->header.time; }
    
    //@{
    //! Particle access
    virtual const Particle* firstParticle ();
    virtual const Particle* nextParticle();
    //@}
  };
#endif  

  typedef std::shared_ptr<ParticleReader> PRptr;
}
  
#endif
  
  
  
