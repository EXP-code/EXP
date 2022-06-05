#include <expand.H>
#include <global.H>
#include <localmpi.H>
#include <ExternalCollection.H>

class UserTest : public ExternalForce
{
private:
  
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();
  int instance;
  static int total;

public:

  UserTest(const YAML::Node& conf);
  ~UserTest();

};

int UserTest::total = 0;

UserTest::UserTest(const YAML::Node& conf) : ExternalForce(conf)
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
  ExternalForce *makerTest(const YAML::Node& conf)
  {
    return new UserTest(conf);
  }
}

class proxytest { 
public:
  proxytest()
  {
    factory["usertest"] = makerTest;
  }
};

static proxytest p;
