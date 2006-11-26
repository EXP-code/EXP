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

  int begin, end;
  vector<int> send, recv, ncount, last;

  void initialize();
  void userinfo();

  MPI_Status status;
  
  void compute_list();
				// For debugging
  void check_body(int);
  void check_send();
  void check_recv();
  void print_send();
  void print_recv();

public:

  //! Constructor
  SatFix(string &line);

  //! Destructor
  ~SatFix();

  //! Get force
  void get_acceleration_and_potential(Component* C);

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
