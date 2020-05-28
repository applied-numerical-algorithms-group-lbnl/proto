
#include <stdlib.h>
#include <stdio.h>
#include <functional>

__global__
void init(int n, int* a)
{
  for(int idx = 0; idx < n; idx +=1)
    {
      //      a[idx]=idx;
      a[idx]=4;
    }

}


int main(int argc, char** argv) 
{
  int n = 16;
  int* aye;
  protoMallocManaged(&aye, n*sizeof(int));

  init<<<1, 1>>>(n, aye);

  printf("out of init\n");

  protoDeviceSynchronize();
  int a0 = aye[0];
  printf("a0=%d\n", a0);
  for(int i=0; i<n; ++i) 
    {

      printf("i = %d, a= %d \n",i,  aye[i]);
    }

  protoFree(aye);

  return 0;
}
