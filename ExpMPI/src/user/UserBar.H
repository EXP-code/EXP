#ifndef _UserBar_H
#define _UserBar_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Bar rotating at specified fraction of corotation

    @param length is the bar length
    @param amp is the bar amplitude 
    @param Ton is the time at the center of the <code>erf</code> turn on
    @param DeltaT is the speard of the turn on
    @param Fcorot is fraction of the length relative to the corotation radius
    @param fixed if nonzero, pattern speed is computed at every step from the current potential, if zero, pattern speed is evaluated once to start.
*/
class UserBar : public ExternalForce
{
private:
  
  string com_name;
  Component *c0;

  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double length, amplitude, Ton, DeltaT, Fcorot;
  bool fixed;

  void userinfo();

public:

  UserBar(string &line);
  ~UserBar();

};

#endif
