#ifndef _ScatterMFP_H
#define _ScatterMFP_H

#include <utility>
#include <random>

#include <ExternalForce.H>

typedef pair<double, int> rpair;

bool less_rpair(const rpair& one, const rpair& two);

//! Compute acceleration based on mean free path (dark matter pressure)
/*! Will probably never use this again . . . */
class ScatterMFP : public ExternalForce
{

private:

  double tau, tauscat, rmax;
  int tautab, nscat;
  double *dtau, *dtau1;
  double dr;
  vector<int> cntr;

  vector <rpair> rr2;

  //@{
  //! Pseudorandom generation
  std::uniform_real_distribution<> unif;
  std::normal_distribution<>       norm;
  //@}

  int cntacc;

  void initialize();
  Component *c;			// Component ptr
  string comp_id;		// Component id
  int mfp_index;		// Index in double attribute list

  // The usual required members

  void get_acceleration_and_potential(Component* C);
  void * determine_acceleration_and_potential_thread(void * arg);

public:

  //! Constructor
  ScatterMFP(const YAML::Node& conf);

  //! Destructor
  ~ScatterMFP();

};


#endif
