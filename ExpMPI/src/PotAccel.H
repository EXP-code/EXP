
#ifndef PotAccel_H
#define PotAccel_H

#include <pthread.h>
#include <string>
#include <list>

#include <Particle.H>
#include <StringTok.H>

class PotAccel;
class Component;

using namespace std;

//! Helper class to pass info to threaded member
struct thrd_pass_PotAccel {
  //! Pointer to force
  PotAccel *t;
  //! Call coefficient thread member (true) or acceleration/potential (false)
  bool coef;
  //! Thread counter id
  int id;
};

//! Defines the basic components of the Force class
/*!
  This class provides the threading members.  Derived classes
  must provide the per thread members for evaluating coeffients
  (if any) and acceleration. 
*/
class PotAccel
{

private:

  // Threading stuff
  thrd_pass_PotAccel *td;
  pthread_t *t;

protected:

  //! Contains parameter database
  list< pair<string, string> > namevalue;
  /*! Call to return <code>value</code> for requested parameter 
    <code>name</code> */
  int get_value(const string& name, string& value);

  //! Pointer to phase space
  vector<Particle> *particles;

  //! Number of bodies
  int nbodies;
  
  //! "Back" pointer to component that instantiated this force
  Component *component;

  //! Used by derived class to initialize any storage and parameters
  virtual void initialize(void) = 0;

  //! Number of particles used on last coefficient evaluation
  int used;

  //! Set flag to add to external potential field
  int use_external;

  //! Coefficient dump flag (true if a useful dump_coefs member is defined)
  bool coef_dump;

public:

  //! Id string
  string id;

  //! Possible geometries
  enum Geometry {sphere, cylinder, cube, slab, other};

  //! sphere, cylinder, cube, slab, other
  Geometry geometry;

  //! Constructor
  /*! <param> line is parsed and info is made available to derived 
    classes on  request */
  PotAccel(string& line);

  //! Destructor
  virtual ~PotAccel();

  //! The main force evaluation routine to be called by clients
  virtual void get_acceleration_and_potential(vector<Particle>*) = 0;

  /*! Used by get_acceleration_and_potential to request that the
    expansion be performed. */
  virtual void determine_coefficients(void) = 0;

  //! Multithreading implementation of the expansion computation
  virtual void * determine_coefficients_thread(void * arg) = 0;

  //! Multithreading implementation of the force computation
  virtual void * determine_acceleration_and_potential_thread(void * arg) = 0;

  /*! Mutex used to be used in the threading implementations to lock
    PCA variance matrix for updating */
  static pthread_mutex_t cc_lock;

  //! Flag to indicate that PCA variance matrix should be recomputed
  static int compute;

  //! For access to parent component
  virtual void RegisterComponent(Component *c) { component = c; }

  //! Thread particle counter
  int *use;

  //! Make some threads
  void exp_thread_fork(bool coef);

  //! Make a mutex
  void make_mutex(pthread_mutex_t *m, char *caller, char *name);

  //! Kill a mutex
  void kill_mutex(pthread_mutex_t *m, char *caller, char *name);

  //! Error
  void bomb(const string& msg);

  //! Return umber of particles used on last coefficient evaluation
  int Used() { return used; }

  //! Set external flag for current force call 
  void SetExternal() { use_external = true; }

  //! Clear external flag
  void ClearExternal() { use_external = false; }

  //! Dump coefficients for this force?
  bool HaveCoefDump() { return coef_dump; }

  //! Dump current coefficients
  virtual void dump_coefs(ostream &out) {};

};

#endif
