#include <math.h>
#include "SatelliteOrbit.h"
#include "expand.h"

/**
   Force dual satellites to be inversion symmetric about the origin
   and to follow a particular analytic orbit
*/
class SatFixOrb : public ExternalForce
{
private:
  
  // Parameters
  bool verbose;
  bool debug;
  double toffset;
  string comp_nam, config;
  int tag;

  Component *c0;

  int begin, end;
  vector<int> send, recv, ncount, last;

  SatelliteOrbit *orb;

  void initialize();
  void userinfo();

  MPI_Status status;
  
  void compute_list();
  void enforce();
				// For debugging
  void check_body(int);
  void check_send();
  void check_recv();
  void print_send();
  void print_recv();

public:

  //! Constructor
  SatFixOrb(string &line);

  //! Destructor
  ~SatFixOrb();

  //! Get force
  void get_acceleration_and_potential(Component* C);

  //! Required threaded member
  void * determine_acceleration_and_potential_thread(void * arg) {return NULL;}
};


extern "C" {
  ExternalForce *makerSatFixOrb(string& line)
  {
    return new SatFixOrb(line);
  }
}

//! Register factory
class proxysatorb { 
public:
  //! Constructor
  proxysatorb()
  {
    factory["satfixorb"] = makerSatFixOrb;
  }
};

proxysatorb p;
