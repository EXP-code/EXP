#ifndef _tidalField_H
#define _tidalField_H

//! Compute tidal acceleration and potential in the Hills limit.
/*!
  The calculation is done in a frame having its cartesian axis 
  directions fixed in space; i.e., the system is *not* phase
  locked with the perturber as is usually the case.
  Unfortunately, I know of no conserved quantities in this frame.
  <br>
  <tt>int freeze_particle(int i)</tt>
  <br>
  Freeze the i-th particle if it has gone too far away from the origin.
*/
class tidalField : public ExternalForce
{
private:
  
  void initialize();
  double hills_omega, hills_p;

  void * determine_acceleration_and_potential_thread(void * arg);

public:

  //! Constructor
  tidalField(string& line);

  //! Destructor
  virtual ~tidalField() {}

};


#endif
