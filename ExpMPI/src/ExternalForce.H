#ifndef _ExternalForce_H
#define _ExternalForce_H

#include <PotAccel.H>

//! External force PotAccel members
/*! 
  This does little but make sure that user is not loading a class
  that is not designed to be an external force.
 */
class ExternalForce : public PotAccel
{
private:

  void determine_coefficients(void);
  void determine_acceleration_and_potential(void);

  void * determine_coefficients_thread(void * arg);

public:
  //! Name of external force (mnemonic)
  string id;
  
  //! Constructor
  /** <param> line passes any parameters to the method */
  ExternalForce(string& line);

  //! Destructor
  virtual ~ExternalForce() {}

  //! Will be called by clients to evaluate force
  void get_acceleration_and_potential(vector<Particle>*);

  /** Multithreaded method (must be thread safe) for force computation
      that must be provided by all derived classes (hence pure virtual) */
  virtual void * determine_acceleration_and_potential_thread(void * arg) = 0;

};

typedef ExternalForce *maker_t(string&);
// extern map<string, maker_t *, less<string> > factory;

#endif

