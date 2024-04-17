#ifndef ParticleFerry_H
#define ParticleFerry_H

#include <vector>
#include <memory>

#include "localmpi.H"
#include "Particle.H"

using namespace std;

//! For interprocess communication, need static structure for MPI mapping
//! 
class ParticleFerry
{
private:

  int nimax, ndmax;
  size_t bufsiz;
  std::vector<char> buf;

  unsigned bufpos, ibufcount, itotcount;
  unsigned _to, _from, _total;

  int keypos, treepos, idxpos;

  void BufferSend();
  void BufferRecv();

  //! Determine size of buffer needed
  void particleBufInit();

public:

  //! Constructor
  ParticleFerry(int nimax, int ndmax);

  //! Destructor
  ~ParticleFerry();

  //! Pack particle into buffer.  Caller provides buffer.
  void particlePack(PartPtr in, char* buffer);

  //! Unpack buffer into particle.  Caller provides buffer.
  void particleUnpack(PartPtr out, char* buffer);

  //! Send message to receiver: get ready for bulk shipment of particles
  void ShipParticles(unsigned to, unsigned from, unsigned& total);

  //@{
  //! Send and receive particles.  Uses ParticleFerry internally
  //! maintained buffer

  void SendParticle(PartPtr part);

  PartPtr RecvParticle();
  //@}

  //! Size needed for a single particle
  size_t getBufsize() { return bufsiz; }
};

typedef std::shared_ptr<ParticleFerry> ParticleFerryPtr;

#endif
