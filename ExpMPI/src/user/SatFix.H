#include <math.h>
#include "expand.h"

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

  SatFix(string &line);
  ~SatFix();

  void get_acceleration_and_potential(vector<Particle>* P);
  void * determine_acceleration_and_potential_thread(void * arg) {return NULL;}
};


extern "C" {
  ExternalForce *makerSatFix(string& line)
  {
    return new SatFix(line);
  }
}

class proxysat { 
public:
  proxysat()
  {
    factory["satfix"] = makerSatFix;
  }
};

proxysat p;
