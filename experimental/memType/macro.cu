/*#define GPU(name,arg) __global__ \
void gpu(name)

#define CPU(name,arg) void cpui##(name)

#define function
*/

__global__
void gpu(double* a_in, unsigned int a_size)
{
  unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
  if(id<a_size)
   a_in[id] += 1;
}

void cpu(double* a_in, unsigned int a_size)
{
  for(int id = 0 ; id < a_size ; id++)
    a_in[id] += 1;
}

template<bool U> 
struct add{
};

template<>
struct add<true>
{
  template<typename... T>
  void operator()(unsigned int nbBlocks, unsigned int nbThreads, T... args)
  {
    gpu<<<nbBlocks,nbThreads>>>(args...);
  }
};
add<true> gpuadd;

template<>
struct add<false>
{
  template<typename... T>
  void operator()(unsigned int nbBlocks, unsigned int nbThreads, T... args)
  {
    cpu(args...);
  }
};
add<false> cpuadd;


#define protoLaunchKernel(memType, Ker, nbBlocks, nbThreads, args...) if(memType) gpu##Ker(nbBlocks,nbThreads,args);\
						else cpu##Ker(nbBlocks,nbThreads,args);

#define protoMalloc(memType,ptr,size) if(memType) {cudaMalloc(&ptr,size);} \
			      else { ptr = (decltype(ptr)) malloc(size); } 
				      


int main()
{
  double* host;
  double* device;

  const bool gpu = true;
  const bool cpu = false;

  unsigned int size = 16;
  unsigned int nBytes = size*sizeof(double);

  protoMalloc(cpu, host, nBytes); 
  protoMalloc(gpu, device,nBytes); 
  
  // just to try, no init

  protoLaunchKernel(cpu,add,1,size,host,size);
  protoLaunchKernel(gpu,add,1,size,device,size);

  return 0;
}
