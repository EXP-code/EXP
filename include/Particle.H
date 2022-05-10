#ifndef Particle_H
#define Particle_H

#include <unordered_map>
#include <vector>
#include <memory>

using namespace std;

// Helper class for buffered binary writes
class ParticleBuffer;

//! Keeps track of all info for one particle
/*!
  The iattrib and dattrib vectors are used by individual components to
  carry additional parameters specific to different particle types.
 */
class Particle
{
private:

  //! MPI buffer
  std::vector<char> buffer;

public:

  //! Default effort value
  static float effort_default;

  //! Mass of particle
  double mass;

  //! Position vector
  double pos[3];

  //! Velocity vector
  double vel[3];

  //! Acceleration vector
  double acc[3];

  //! Gravitational potential due to particles  
  double pot;

  //! Gravitational potential from external forces
  double potext;
  
  //! Integer attributes
  vector<int> iattrib;

  //! Real (double) attributes
  vector<double> dattrib;

  //! Multistep level
  unsigned level;

  //! Internal process time step (<0 means ignore)
  float dtreq;

  //! Internal size scale (<0 means ignore)
  float scale;

  //! Scale factor for partitioning work
  float effort;

  //! Index for ordering
  unsigned long indx;

  //! Tree key
  unsigned tree;

  //! Hash key
  unsigned long key;
  
  typedef std::pair<unsigned short, unsigned short> speciesKey;
  static const speciesKey defaultKey;

  //! Species key (this is a cache value to prevent endless on-the-fly
  //! reconstruction of the species key)
  speciesKey skey;

  //! Constructor
  Particle();

  //! Constructor with presized attribute lists
  Particle(unsigned niatr, unsigned ndatr);

  //! Copy constructor
  Particle(const Particle &);

  //! Read particles from file
  void readAscii(bool indexing, int seq, std::istream* fin);

  //! Read particles from file 
  void readBinary(unsigned rsize, bool indexing, int seq, std::istream *in);

  //! Write a particle in ascii format
  void writeAscii(bool indexing, bool accel, std::ostream* out);

  //! Write a particle in binary format (PSP)
  void writeBinaryBuffered(unsigned rsize, bool indexing, std::ostream *out, ParticleBuffer& buf) const;

  //! Write a particle in binary format (PSP), threaded version
  void writeBinaryThreaded
  (unsigned rsize, bool indexing, ParticleBuffer& buf, int count) const;

  //! Write a particle in binary format (PSP)
  void writeBinary(unsigned rsize, bool indexing, std::ostream *out) const;

  //! Write a particle in binary format (PSP) to buffer for MPI
  int writeBinaryMPI(char* buf, unsigned rsize, bool indexing);
  
  //! Particle buffer size
  unsigned getMPIBufSize(unsigned rsize, bool indexing)
  {
    unsigned csize = (8 + dattrib.size()) * rsize + iattrib.size() * sizeof(int);
    if (indexing) csize += sizeof(unsigned long);
    return csize;
  }
};

//! Helper class for writing a buffer of particles
class ParticleBuffer
{
public:

  //! Buffer maximum particle count
  const unsigned maxBufCount = 32768;

protected:

  //! Write buffer
  std::vector<char> charBuffer;

  //! Particle size
  unsigned parSize;

  //! Buffer size
  unsigned bufSize;

  //! Buffer particle count
  unsigned bufCount;

public:

  //! Buffer location
  char* Loc;

  //! Constructor
  ParticleBuffer(unsigned rsize, bool indexing, const Particle *p);

  //! Write buffer
  void writeBuffer(std::ostream *out, bool finish=false);

  //! Increment the particle counter
  void operator++(int) { bufCount++; }

  //! Raw access to particle buffer
  char* operator()() { return charBuffer.data(); }

  //! Set location to the nth particle
  void setLocation(int n) { Loc = charBuffer.data() + parSize*n; }

  //! Get memory per particle
  unsigned pMemSize() { return parSize; }

};


typedef std::shared_ptr<Particle> PartPtr;
typedef std::unordered_map<unsigned long, PartPtr> PartMap;
typedef std::unordered_map<unsigned long, PartPtr>::iterator PartMapItr;
typedef std::pair<unsigned long, PartPtr> PMapType;
std::ostream& operator<< (std::ostream& os, const PMapType& p);

#endif
