#ifndef ParticleFerry_H
#define ParticleFerry_H

#include <vector>
#include "localmpi.h"
#include "Particle.H"

using namespace std;

//! Maximum number of *extra* integer components
const static int nimax = 5;

//! Maximum number of *extra* double components
const static int ndmax = 16;

//! For interprocess communication, need static structure for MPI mapping
class Partstruct
{
public:
  /// Constructor
  Partstruct();

  ///< Type:                  Offsets:
  double mass;                ///< 0        0
  double pos[3];              ///< 1        1,2,3
  double vel[3];              ///< 2        4,5,6
  double acc[3];              ///< 3        7,8.9
  double pot;                 ///< 4        10
  double potext;              ///< 5        11
  float  dtreq;               ///< 6        12
  float  scale;               ///< 7        13
  float  effort;              ///< 8        14
  unsigned level;             ///< 9        15
  unsigned long indx;         ///< 10       16
  unsigned tree;	      ///< 11       17
#ifdef I128
  uint128 key;		      ///< 12       18
#else
  unsigned long key;	      ///< 12       18
#endif
  unsigned nicnt;             ///< 13       19
  unsigned ndcnt;             ///< 14       20 (e.g. for nimax=5 and ndmax=12)
  int    iatr[nimax];         ///< 15       21,22,23,24,25
  double datr[ndmax];         ///< 16       26,27,28,29,30,31,32,33,34,35,37,38
};

class ParticleFerry
{
private:

  Partstruct *buf;
  unsigned ibufcount, itotcount;
  unsigned _to, _from, _total;
  unsigned pk_lo, pk_hi;
#ifdef I128
  uint128  key_lo, key_hi;
#else
  unsigned long key_lo, key_hi;
#endif

  void BufferSend();
  void BufferRecv();
  void bufferKeyCheck();	// For debugging
  void PartstructInit(Partstruct*);

public:

  ParticleFerry();

  ~ParticleFerry();

  void ShipParticles(unsigned to, unsigned from, unsigned& total);

#ifdef I128
  void SendParticle(Particle& part, unsigned seq,  uint128  key);

  bool RecvParticle(Particle& part, unsigned& seq, uint128& key);
#else
  void SendParticle(Particle& part, unsigned seq,  unsigned long  key);

  bool RecvParticle(Particle& part, unsigned& seq, unsigned long& key);
#endif

  void SendParticle(Particle& part);

  bool RecvParticle(Particle& part);
  
  void part_to_Particle(Partstruct& str, Particle& cls);

  void Particle_to_part(Partstruct& str, Particle& cls);

  MPI_Datatype Particletype;
};

#endif
