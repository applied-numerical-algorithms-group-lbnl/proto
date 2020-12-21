#define GPU(name) __global__ \
void gpu_##name

#define CPU(name) void cpu_##name

#define FUNCTOR_BEGIN(name) template<bool U> struct base##name{\
};\
\
template<> struct base##name<true>\
{\
  template<typename... T>\
  void operator()(unsigned int nbBlocks, unsigned int nbThreads, T... args)\
  {\
    gpu_##name<<<nbBlocks,nbThreads>>>(args...);\
  }\
};\
template<> struct base##name<false>\
{\
  template<typename... T>\
  void operator()(unsigned int nbBlocks, unsigned int nbThreads, T... args)\
  {\
    cpu_##name(args...);\
  }\
};\
base##name<true> gpu##name;\
base##name<false> cpu##name;


GPU(add)(double* a_in, unsigned int a_size)
{
  unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
  if(id<a_size)
   a_in[id] += 1;
};

CPU(add)(double* a_in, unsigned int a_size)
{
  for(int id = 0 ; id < a_size ; id++)
    a_in[id] += 1;
};

FUNCTOR(add);


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
  
  protoLaunchKernel(cpu,add,1,size,host,size);
  protoLaunchKernel(gpu,add,1,size,device,size);

  return 0;
}
