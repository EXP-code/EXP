#include <expand.h>
#include <global.H>
#include <localmpi.h>
#include <ExternalCollection.H>

class UserTest : public ExternalForce
{
private:
  
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();
  int instance;
  static int total;

public:

  UserTest(string &line);
  ~UserTest();

};

int UserTest::total = 0;

UserTest::UserTest(string &line) : ExternalForce(line)
{
  id = "Test";
  total++;
  instance = total;
  if (myid==0) 
    cout << "Just made a UserTest! Instance=" << instance << endl;
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
  cout << "Process " << myid 
       << ", id=" << id << ": Time=" << tnow 
       << ", Instance=" << instance << endl;
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
