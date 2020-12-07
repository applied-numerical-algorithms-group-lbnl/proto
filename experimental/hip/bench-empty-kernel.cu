#include<iostream>
#include<Proto_gpu.H>

__global__
void add(double *ptr1, double *ptr2, unsigned int n)
{
  unsigned int id = threadIdx.x + blockIdx.x * blockDim.x;
  if(id<n)
  {
    ptr1[id] += ptr2[id];
  }
}


__global__
void copy(double *ptr1, double *ptr2, unsigned int n)
{
  unsigned int id = threadIdx.x + blockIdx.x * blockDim.x;
  if(id<n)
  {
    ptr1[id] = ptr2[id];
  }
}

__global__
void init(double *ptr1, double *ptr2, unsigned int n)
{
  unsigned int id = threadIdx.x + blockIdx.x * blockDim.x;
  if(id<n)
  {
    ptr1[id] = 1;
    ptr2[id] = 2;
  } 
}

int main()
{
//  cudaSetDevice(3);
  double *p1;
  double *p2;

  unsigned int n = std::pow(10,7);
  protoDeviceSynchronize();

  protoMalloc(p1, n*sizeof(double));
  protoMalloc(p2, n*sizeof(double));

  init<<<n+256-1,256>>>(p1,p2,n); 

  double results[5*2];
  protoEvent_t start[5*2], stop[5*2];

  for(int i = 0 ; i < 10 ; i++)
  {	
    protoEventCreate(&start[i]);
    protoEventCreate(&stop[i]);
  }

  int acc = 0;
  for(int size = 10 ; size <= n ; size *=10)
  {
    unsigned int nbThread = 512;
    unsigned int nbBlock  = (size + nbThread - 1) / nbThread;
    protoEventRecord(start[acc]);
    for(int i = 0 ; i < 1000 ; i++)
    {
      copy<<<nbBlock,nbThread>>>(p1,p2,size);
    }
    protoDeviceSynchronize();
    protoEventRecord(stop[acc]);
    acc++;
  }
  std::cout << acc << std::endl;

  for(int i = 0 ; i < acc ; i++)
  {
    float milliseconds = 0;
    protoDeviceSynchronize();
    protoEventElapsedTime(&milliseconds, start[i], stop[i]);
    results[i] = milliseconds;
  }

  std::cout << " \\addplot plot coordinates {" <<std::endl;
  
  unsigned int j = 10;
  for(int i = 0 ; i < acc ; i++)
  {
    std::cout << "( " <<j << " , " << results[i] << ") \n";
    j *= 10;
  }
  std::cout << "};" <<  std::endl;

  return 0;
}
