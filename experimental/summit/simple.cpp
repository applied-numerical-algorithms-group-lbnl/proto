
#include <stdlib.h>
#include <stdio.h>
#include <functional>

void init(int n, int* a)
{
  for(int idx = 0; idx < n; idx +=1)
    {
      a[idx]=idx;
    }

}


int main(int argc, char** argv) 
{
  int n = 16;
  int* aye;
  aye = new int[n];

  init(n, aye);

  printf("out of init\n");

  int a0 = aye[0];
  printf("a0=%d\n", a0);
  for(int i=0; i<n; ++i) 
    {
      printf("i = %d, a= %d \n",i,  aye[i]);
    }

  delete[] aye;

  return 0;
}
