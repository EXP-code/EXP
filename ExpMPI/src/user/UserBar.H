#ifndef _UserBar_H
#define _UserBar_H

#include <Particle.H>
#include <AxisymmetricBasis.H>
#include <ExternalCollection.H>

/** Bar rotating at specified fraction of corotation

    @param length is the bar length
    @param amp is the bar amplitude 
    @param Ton is the time at the center of the <code>erf</code> turn on
    @param Toff is the time at the center of the <code>erf</code> turn off
    @param DeltaT is the spread of the turn on
    @param Fcorot is fraction of the length relative to the corotation radius
    @param fixed if nonzero, pattern speed is computed at every step from the current potential, if zero, pattern speed is evaluated once to start (default: false).
    @param soft use the shallow form of the potential (default: false)
*/
class UserBar : public ExternalForce
{
private:
  
  string com_name;
  Component *c0;

  void determine_acceleration_and_potential(void);
  void determine_acceleration_and_potential(vector<Particle>*);
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

  double length, amplitude, Ton, Toff, DeltaT, Fcorot;
  bool fixed, soft;
  string filename;

  static const double numfac = 3.86274202023190e-01;
  bool firstime;
  double posang, lastomega, lasttime;
  double Lz, Lz0, Lzbar, omega, Iz, afac, b5;
  string name;

  void userinfo();

public:

  //! Constructor
  UserBar(string &line);

  //! Destructor
  ~UserBar();

};

#endif
