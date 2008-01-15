#ifndef ParticleFerry_H
#define ParticleFerry_H

#include <vector>
#include "localmpi.h"
#include "Particle.H"

using namespace std;

//! Maximum number of *extra* integer components
const static int nimax = 5;
//! Maximum number of *extra* double components
const static int ndmax = 12;

//! For interprocess communication, need static structure for MPI mapping
struct Partstruct
{
  ///< Type:                  Offsets:
  double mass;                ///< 0        0
  double pos[3];              ///< 1        1,2,3
  double vel[3];              ///< 2        4,5,6
  double acc[3];              ///< 3        7,8.9
  double pot;                 ///< 4        10
  double potext;              ///< 5        11
  float  dtreq;               ///< 6        12
  unsigned level;             ///< 7        13
  unsigned long indx;         ///< 8        14
  unsigned long long key;     ///< 9        15
  unsigned nicnt;             ///< 10       16
  unsigned ndcnt;             ///< 11       17
  int    iatr[nimax];         ///< 12       18,19,20,21,22 (e.g.)
  double datr[ndmax];         ///< 13       23,24,25,26,27 (e.g.)
};

class ParticleFerry
{
private:

  Partstruct *buf;
  unsigned ibufcount, itotcount;
  unsigned _to, _from, _total;

  void BufferSend();
  void BufferRecv();
  void bufferKeyCheck();	// For debugging

public:

  static MPI_Datatype Particletype;

  ParticleFerry();

  void ShipParticles(int to, int from, unsigned& total);

  void SendParticle(Particle& part, unsigned seq, unsigned long long key);

  bool RecvParticle(Particle& part, unsigned& seq, unsigned long long& key);

  void SendParticle(Particle& part);

  bool RecvParticle(Particle& part);
  
  void part_to_Particle(Partstruct& str, Particle& cls);

  void Particle_to_part(Partstruct& str, Particle& cls);
};

#endif
