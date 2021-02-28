#include <cudaUtil.cuH>

#include "UserTestCuda.H"


__global__ void cuda_hello(int myid, int id, double time, int count)
{
  printf("GPU says: Process %d, id=%d: Time=%f, Instance=%d\n", myid, id, time, count);
}

void UserTestCuda::cuda_user_test(int myid, int id, double time, int count)
{
  cuda_hello<<<1,1>>>(myid, id, time, count); 
}
