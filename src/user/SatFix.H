#include <cmath>
#include <expand.H>
#include <ExternalCollection.H>

//! Force dual satellites to be inversion symmetric about the origin
class SatFix : public ExternalForce
{
private:
  
  // Parameters
  bool verbose;
  bool debug;
  string comp_name;

  Component *c0;

  vector<int> owner;

  void initialize();
  void userinfo();

  MPI_Status status;
  
  void compute_list();
				// For debugging
  void check_body(int);
  void print_body(int, string&);

  //! Valid keys for YAML configurations
  static const std::set<std::string> valid_keys;

public:

  //! Constructor
  SatFix(const YAML::Node& conf);

  //! Destructor
  ~SatFix();

  //! Get force
  void get_acceleration_and_potential(Component* C);

  //! Required threaded member
  void * determine_acceleration_and_potential_thread(void * arg) {return NULL;}
};


extern "C" {
  ExternalForce *makerSatFix(const YAML::Node& conf)
  {
    return new SatFix(conf);
  }
}

//! Register factory
class proxysatfix { 
public:
  //! Constructor
  proxysatfix()
  {
    factory["satfix"] = makerSatFix;
  }
};

proxysatfix p;
