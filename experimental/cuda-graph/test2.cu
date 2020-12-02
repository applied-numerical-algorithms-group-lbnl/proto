#include <iostream>

#define N 500000 
#define NSTEP 1000
#define NKERNEL 20

__global__ 
void shortKernel(float * out_d, float * in_d){
  int idx=blockIdx.x*blockDim.x+threadIdx.x;
  if(idx<N) out_d[idx]=1.23*in_d[idx];
}

__global__
void bigKernel(float * out_d, float * in_d){
  int idx=blockIdx.x*blockDim.x+threadIdx.x;
  if(idx<N)
  {
    for(int it = 0; it<NKERNEL ; it++)
      out_d[idx]=1.23*in_d[idx];
  }
}

__global__ 
void initKernel(float * out_d, float * in_d){
  int idx=blockIdx.x*blockDim.x+threadIdx.x;
  if(idx<N) out_d[idx]=0;
  if(idx<N) in_d[idx] =1;
}

int main()
{
  float * in_d, *out_d;
  cudaMalloc(&in_d, N*sizeof(float));
  cudaMalloc(&out_d, N*sizeof(float));

  unsigned int threads = 256;
  unsigned int blocks = (N+threads-1)/threads;

  initKernel<<<blocks, threads, 0, 0>>>(out_d, in_d);

  std::cout << " classic " << std::endl;
  cudaEvent_t start, stop;

  {
    float milliseconds = 0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    // start CPU wallclock timer
    for(int istep=0; istep<NSTEP; istep++){
      for(int ikrnl=0; ikrnl<NKERNEL; ikrnl++){
        shortKernel<<<blocks, threads, 0, 0>>>(out_d, in_d);
        cudaStreamSynchronize(0);
      }
    }
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << " time: " << milliseconds << " ms" << std::endl;
  }
  //end CPU wallclock time

  std::cout << " graph " << std::endl;
  initKernel<<<blocks, threads, 0, 0>>>(out_d, in_d);
  bool graphCreated=false;
  cudaGraph_t graph;
  cudaGraphExec_t instance;
  {
    cudaStream_t stream;
    cudaStreamCreate ( &stream);
    float milliseconds = 0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    for(int istep=0; istep<NSTEP; istep++){
      if(!graphCreated)
      {
        cudaStreamBeginCapture(stream, cudaStreamCaptureModeGlobal);
        for(int ikrnl=0; ikrnl<NKERNEL; ikrnl++){
          shortKernel<<<blocks, threads, 0, stream>>>(out_d, in_d);
        }
        cudaStreamEndCapture(stream, &graph);
        cudaGraphInstantiate(&instance, graph, NULL, NULL, 0);
        graphCreated=true;
      }
      cudaGraphLaunch(instance, stream);
      cudaStreamSynchronize(stream);
    }
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << " time: " << milliseconds << " ms" << std::endl;
  }

  std::cout << " Fusion " << std::endl;
  initKernel<<<blocks, threads, 0, 0>>>(out_d, in_d);
  {
    float milliseconds = 0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);

    // start CPU wallclock timer
    for(int istep=0; istep<NSTEP; istep++){
      bigKernel<<<blocks, threads, 0, 0>>>(out_d, in_d);
      cudaStreamSynchronize(0);
    }
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&milliseconds, start, stop);
    std::cout << " time: " << milliseconds << " ms" << std::endl;
  }

  return 0;
}
