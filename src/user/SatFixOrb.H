#include <memory>
#include <cmath>

#include <expand.H>
#include <SatelliteOrbit.h>

#include <ExternalCollection.H>

/**
   Force dual satellites to be inversion symmetric about the origin
   and to follow a particular analytic orbit
   
   The particle with integer attribute <code>tag</code> will be forced
   to maintain the analytic orbit defined by SatelliteOrbit with
   parameters specified in the file <code>config</code>.

   The remaining particles will be forced to be inversion symmetric by
   adjacent pairs.
   

   @param tag is the integer attribute of the particle whose orbit
   will be enforced

   @param config has the parameters for SatelliteOrbit

   @param compname is the point mass component containing the
   particles to be tracked

   @param toffset shifts the origin of the orbit at T=0 from the
   orbit's pericenter
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
  std::vector<unsigned int> ncount, last;
  std::vector<int> send, recv;

  std::shared_ptr<SatelliteOrbit> orb;

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

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  SatFixOrb(const YAML::Node& conf);

  //! Destructor
  ~SatFixOrb();

  //! Get force
  void get_acceleration_and_potential(Component* C);

  //! Required threaded member
  void * determine_acceleration_and_potential_thread(void * arg) {return NULL;}
};


extern "C" {
  ExternalForce *makerSatFixOrb(const YAML::Node& conf)
  {
    return new SatFixOrb(conf);
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
