#ifndef _generateRelaxation_H
#define _generateRelaxation_H

//! Compute relaxation diagnoistics
class generateRelaxation : public ExternalForce
{

private:
  
  int epos;
  int done;
  void initialize();

  void * determine_acceleration_and_potential_thread(void * arg);

public:

  //! Constructor
  generateRelaxation(const YAML::Node& conf);

};

#endif
