#ifndef _ExternalCollection_H
#define _ExternalCollection_H

#include <map>

#include "ComponentContainer.H"
#include "ExternalForce.H"

using namespace std;

//! factory for dynamically loaded forces
extern map<string, maker_t *, less<string> > factory;
extern map<string, maker_t *, less<string> >::iterator fitr;


//! Container class for external force types
/*!
  These are supplied as a file with one type to a line.  The constructor
  does "hand coded" reflection to make a chain of these.  Therefore, the
  external forces will be applied in order of their listing (in case this
  matters).
  <br>
  Format: single line per force, parameters separated from name by ":"
*/
class ExternalCollection
{
private:

  list<ExternalForce *>::iterator sitr;

  list<void *> dl_list;		// list to hold handles for dynamic libs 
  list<void *>::iterator itr; 
  void dynamicload(void);
  vector<string> getlibs(void);

public:
  
  //! List of forces we create
  list<ExternalForce *> force_list;	

  //! Constructor
  ExternalCollection();

  //! Destructor
  ~ExternalCollection();
  
  //! Reads configuration database and initializes requested external forces
  void initialize();

  //! Provides the same interface to compute acceleration as ExternalForce
  void get_acceleration_and_potential(vector<Particle>*);

  //! Call the finish member for all external forces
  void finish();
};

#endif
