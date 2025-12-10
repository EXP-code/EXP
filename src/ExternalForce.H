#ifndef _ExternalForce_H
#define _ExternalForce_H

#include "PotAccel.H"
#include "Component.H"

//! External force PotAccel members
/*! 
  This does little but make sure that user is not loading a class
  that is not designed to be an external force.
 */
class ExternalForce : public PotAccel
{
private:

  void determine_coefficients(void);

  void * determine_coefficients_thread(void * arg);

protected:

  /** Called to evaluate acceleration
      Default version calls the provided thread function.
      Can be overloaded to provide additional non-threaded
      tasks.
  */
  virtual void determine_acceleration_and_potential(void);

  /** User info divider
      Make demarcation in output log
  */
  void print_divider(void);

public:

  //! Name of external force (mnemonic)
  string id;
  
  //! Constructor
  /** \param line passes any parameters to the method */
  ExternalForce(const YAML::Node& config);

  //! Destructor
  virtual ~ExternalForce() {}

  //! Will be called by clients to evaluate force
  void get_acceleration_and_potential(Component * compIn);

  /** Multithreaded method (must be thread safe) for force computation
      that must be provided by all derived classes (hence pure virtual) */
  virtual void * determine_acceleration_and_potential_thread(void * arg) = 0;

  //! Finish and clean-up (caching data necessary for restart)
  virtual void finish() {}

  // #if HAVE_LIBCUDA==1
#if HAVE_LIBCUDA==1
  //! Copy particles from device for non-cuda forces
  void getParticlesCuda(Component *c);
#endif


};

typedef ExternalForce *maker_t(const YAML::Node&);

#endif

