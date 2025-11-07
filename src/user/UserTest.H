#ifndef _User_Test_H
#define _User_Test_H

#include "ExternalCollection.H"

/** Simple skeleton user routine
    Does nothing but print out messages at each call
*/
class UserTest : public ExternalForce
{
private:
  
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

public:
  //! Constructor
  UserTest(const YAML::Node& conf);

  //! Destructor
  ~UserTest();

};

#endif
