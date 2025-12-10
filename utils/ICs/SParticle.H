#ifndef _SParticle_H
#define _SParticle_H

#include "Particle.H"
#include "localmpi.H"

/**
   An abbreviated particle class for transferring phase space between
   processors and writing initial conditions.  Changes in or additions
   to the Particle class should not affect this class.
*/
class SParticle
{
private:

  //! Copy Particle contents
  void cpy(const Particle& p)
  {
    mass = p.mass;
    for (size_t k=0; k<3; k++) {
      pos[k] = p.pos[k];
      vel[k] = p.vel[k];
    }
  }

public:
  //@{
  //! Data fields
  double mass;
  double pos[3];
  double vel[3];
  //@}

  //! Null constructor
  SParticle() {}

  //! Construct instance from a Particle
  SParticle(const Particle& p) { cpy(p); }

  //! Copy contents of a Particle to an SParticle
  void operator[](const Particle& p) { cpy(p); }

};

/**
   Create an MPI particle type SParticle.  The class creates and
   registers the particle type.  It is freed when the instance goes
   out of scope and/or is destroyed.
*/
class SPtype
{
private:
  MPI_Datatype Particletype;
public:

  SPtype();
  ~SPtype();

  MPI_Datatype operator()() { return Particletype; }
};


#endif
