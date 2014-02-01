#ifndef Particle_H
#define Particle_H

class Particle 
{
 public:

  unsigned level;
  double mass;
  double pos[3];
  double vel[3];
  int Z;
  int C;

  // Constructor
  Particle();

  // Copy constructor
  Particle(const Particle &p);
};

#endif
