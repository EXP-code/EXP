
#include <ExternalCollection.H>

class UserTest : public ExternalForce
{
private:
  
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();

public:

  UserTest(string &line);

};


UserTest::UserTest(string &line) : ExternalForce(line)
{
  cout << "Just made a UserTest!\n";
}

void UserTest::initialize()
{
}

void * determine_acceleration_and_potential_thread(void * arg) 
{
  return NULL;
}


extern "C" {
  ExternalForce *maker(string& line)
  {
    return new UserTest(line);
  }
}

class proxy { 
public:
  proxy()
  {
    factory["usertest"] = maker;
  }
};

proxy p;


