#include <math.h>
#include "expand.h"

//! Force dual satellites to be inversion symmetric about the origin
class SatFix : public ExternalForce
{
private:
  
  // Parameters
  bool verbose;
  string comp_name;

  Component *c0;

  vector<int> plocate;
  int myid0, myid1;

  void initialize();
  void userinfo();

public:

  //! Constructor
  SatFix(string &line);

  //! Destructor
  ~SatFix();

  //! Get force
  void get_acceleration_and_potential(vector<Particle>* P);

  //! Required threaded member
  void * determine_acceleration_and_potential_thread(void * arg) {return NULL;}
};


extern "C" {
  ExternalForce *makerSatFix(string& line)
  {
    return new SatFix(line);
  }
}

//! Register factory
class proxysat { 
public:
  //! Constructor
  proxysat()
  {
    factory["satfix"] = makerSatFix;
  }
};

proxysat p;
