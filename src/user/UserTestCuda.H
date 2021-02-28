#ifndef _UserTestCuda_H
#define _UserTestCuda_H

#include <ExternalForce.H>

class UserTestCuda : public ExternalForce
{
private:
  
  void * determine_acceleration_and_potential_thread(void * arg);
  void initialize();
  int instance, callno, maxcall;
  static int total;

#if HAVE_LIBCUDA == 1
  void cuda_user_test(int myid, int id, double time, int count);
#endif

public:

  UserTestCuda(const YAML::Node& conf);
  ~UserTestCuda();

};

#endif
