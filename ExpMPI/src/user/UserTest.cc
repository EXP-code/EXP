#include "expand.h"

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

UserTest::UserTest(string &line) : ExternalForce(line)
{
  id = "Test";
  cout << "Just made a UserTest!\n";
}

UserTest::~UserTest()
{
}

void UserTest::initialize()
{
}

void * UserTest::determine_acceleration_and_potential_thread(void * arg) 
{
  int id = *((int*)arg);
  cout << "Process " << myid << ", id=" << id << ": Time=" << tnow << endl;
  return (NULL);
}


extern "C" {
  ExternalForce *makerTest(string& line)
  {
    return new UserTest(line);
  }
}

class proxytest { 
public:
  proxytest()
  {
    factory["usertest"] = makerTest;
  }
};

proxytest p;
