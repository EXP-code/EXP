
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
  if (myid==0) cout << "Just made a UserTest! Instance=" << instance << endl;
}

void UserTest::initialize()
{
}

void * determine_acceleration_and_potential_thread(void * arg) 
{
  int id = *((int*)arg);
  if (myid==0 && id==0) cout << "A UserTest[" << instance << "] calls accel!"
			     << endl;
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


