#ifndef _User_Test_H
#define _User_Test_H

#include <ExternalCollection.H>

class UserTest : public ExternalForce
{
private:
  
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

public:

  UserTest(string &line);
  ~UserTest();

};

#endif
