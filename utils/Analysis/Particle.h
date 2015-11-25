#ifndef Particle_H
#define Particle_H

#define _REDUCED

class Particle 
{
 public:

  unsigned level;
  double mass;
  double pos[3];
  double vel[3];

  // Constructor
  Particle();

  // Copy constructor
  Particle(const Particle &p);
};

#endif
