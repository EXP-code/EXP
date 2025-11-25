#include "cudaUtil.cuH"

#include "UserTestCuda.H"


__global__ void cuda_hello(int myid, int id, double time, int count, int num)
{
  printf("GPU says: Process %d, id=%d: Time=%f, Instance=%d, #=%d\n", myid, id, time, count, num);
}

void UserTestCuda::cuda_user_test(int myid, int id, double time, int count, int num)
{
  cuda_hello<<<1,1>>>(myid, id, time, count, num); 
}
