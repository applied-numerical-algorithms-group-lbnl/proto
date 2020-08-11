#include <Proto_gpu.H>
#define DIM2D
__global__
void initOne(double *ptr, unsigned int size)
{
	for(unsigned int i = threadIdx.x +blockDim.x*blockIdx.x ; i < size ; i+= gridDim.x*blockDim.x )
		ptr[i] = 1;
}

__global__
void lapl1DBasic(double *in, double *out, const unsigned int size)
{
	double acc;
	for(unsigned int i = threadIdx.x +blockDim.x*blockIdx.x ; i < size ; i+= blockDim.x * gridDim.x )
	{
		if(i>0 && i < size-1)
		{
			acc = 4*in[i];
			acc-= in[i-1]; 
			acc-= in[i+1]; // 4 for check
			out[i] = acc; //4*in[i] - in[i-1] - in[i+1]; // 4 for checking
		}
	}
}

__global__
void lapl1DBasicV2(double *in, double *out, unsigned int size)
{
	unsigned int i = threadIdx.x +blockDim.x*blockIdx.x; 
	if(i>0 && i < size-1)
	{
		out[i] = 4*in[i] - in[i-1] - in[i+1]; // 4 for check
	}
}

const int NN=256;

__global__
void lapl1DShared(double *in, double *out, const unsigned int size)
{
	//extern __shared__ double mem[];
	__shared__ double mem[NN+2];
	const unsigned int myId = threadIdx.x +1;
	const unsigned int blockSize = gridDim.x*blockDim.x;
	for(unsigned int i = myId + blockDim.x*blockIdx.x ; i < size ; i+= blockSize)
	{
		mem[myId] = in[i];
		if(myId-1==0)
		{
			mem[0] = in[i-1];
			const unsigned int idx = i-1 + blockSize;
			if(idx<size)
			mem[NN+1] = in[idx];
		}

		__syncthreads();

		if(i<size-1)
		out[i] = 4*mem[myId] - mem[myId-1] - mem[myId+1]; // 4 for check
	}
}

__global__
void lapl1DSharedV2(double *in, double *out, const unsigned int size)
{
	//extern __shared__ double mem[];
	__shared__ double mem[NN+2];
	const unsigned int myId = threadIdx.x +1; //myId>0
	const unsigned int blockSize = blockDim.x;
	const unsigned i = myId + blockDim.x*blockIdx.x; // i>1
	if(i<size-1)
	{
	  mem[myId] = in[i];
	  if(myId-1==0)
	  {
			mem[0] = in[i-1];
			const unsigned int idx = i-1 + blockSize;
			if(idx<size)
			mem[NN+1] = in[idx];
	  } 
	}
	double acc;
	__syncthreads();

	acc = 4*mem[myId];
	acc-= mem[myId-1]; 
	acc-= mem[myId+1]; // 4 for check

	if(i<size-1)
		out[i] = acc; //4*in[i] - in[i-1] - in[i+1]; // 4 for checking
}



__global__
void lapl2DCoal(double *in, double *out, const unsigned int size)
{
	//extern __shared__ double mem[];
	//extern __shared__ double mem[]; //[(NN+2)*(NN+2)];
	const unsigned id = threadIdx.x + blockDim.x*blockIdx.x; // i>1

	const int x = id % size;
	const int y = id / size;

	if( x<size - 1 && y < size -1 && x > 0 && y >0)
	{
		double acc = 6*in[id];
		acc-= sin(in[id-1]); 
		acc-= sin(in[id+1]); 
		acc-= exp(in[id + size]); 
		acc-= exp(in[id - size]); 
		out[id] = acc; //4*in[i] - in[i-1] - in[i+1]; // 4 for checking
	}
}

__global__
void lapl2DShared(double *in, double *out, const unsigned int size)
{
	//extern __shared__ double mem[];
	//extern __shared__ double mem[]; //[(NN+2)*(NN+2)];
	__shared__ double mem[18*18];
	const unsigned int myIdx = threadIdx.x ; //myId>0
	const unsigned int myIdy = threadIdx.y ; //myId>0
	const unsigned i = myIdx + blockDim.x*blockIdx.x; // i>1
	const unsigned j = myIdy + blockDim.y*blockIdx.y; // i>1

	const unsigned int id = i + j * size;
	const unsigned int myId = myIdx + 1 + (16+2) * (myIdy+1);

	if(id <size*size)
	{
	  mem[myId] = in[id];
	 /* if(myIdx==0)
	  {
			mem[myId-1] = in[id-1];
			const unsigned int idx = myId + 16; //myId-1+ NN +1
			if(idx<size)
			mem[idx] = in[id + 16];
	  }
           if(myIdy==0)
          {
                        mem[myId-18] = in[id-size];
                        const unsigned int idy = myId + (16)*(16+2); //myId-1+ NN +1
                        if(idy<size)
                        mem[idy] = in[id + size*(16+1)];
          }
 */
	}
	double acc;
	__syncthreads();

	acc = 6*mem[myId];
	acc-= sin(mem[myId-1]); 
	acc-= sin(mem[myId+1]); // 4 for check
	acc-= exp(mem[myId-(16+2)]); // 4 for check
	acc-= exp(mem[myId+16+2]); // 4 for check

	if(i<size - 1 && j < size -1 )
		out[id] = acc; //4*in[i] - in[i-1] - in[i+1]; // 4 for checking
}

__global__
void lapl2D(double *in, double *out, const unsigned int size)
{
	const unsigned i = threadIdx.x + blockDim.x*blockIdx.x; //
	const unsigned j = threadIdx.y + blockDim.y*blockIdx.y; //

	const unsigned int id = i + j * size;

	if(i>0 && j>0 && i<size-1 && j < size-1)
	{
	  double acc;

	  acc = 6*in[id];
	  acc-= sin(in[id-1]); 
	  acc-= sin(in[id+1]); 
	  acc-= exp(in[id-size]);
	  acc-= exp(in[id+size]);
	  
	  out[id] = acc; //4*in[i] - in[i-1] - in[i+1]; // 4 for checking
	}
}


int main()
{

	//int size = 100000000;
	int size = 10000;
	bool good =true;

#ifdef DIM1D
	double * out;
	double * in;
	double * host = new double[size];
	protoMalloc(&in, size * sizeof(double));
	protoMalloc(&out, size * sizeof(double));


	protoLaunchKernel(initOne,1, NN, in,size);	
	protoLaunchKernel(lapl1DBasicV2, (size+NN-1)/NN, NN, in, out, size);


	protoLaunchKernel(initOne,1, NN, in,size);	
	protoLaunchKernel(lapl1DSharedV2, (size+NN-1)/NN, NN, in, out, size);

	protoLaunchKernel(initOne,1, NN, in,size);
	protoLaunchKernel(lapl1DBasic, 256, NN, in, out, size);

	protoLaunchKernel(initOne,1, NN, in,size);	
	protoLaunchKernel(lapl1DShared, 256, NN, in, out, size);

	protoMemcpy(host,out,size*sizeof(double), protoMemcpyDeviceToHost);


	if(host[0] != 1 && host[size-1] != 1) good = false;
	for(int i = 1; i<size-1 ; i++)
		if(host[i] != 2)
		{
			good = false;
			std::cout << " error host["<<i<<"] = " << host[i] << std::endl;
		        break;	
		}
	if(good==true) std::cout << " stencil checked " << std::endl;
#endif
#ifdef DIM2D


	double * out2D;
	double * in2D;
	double * host2D = new double[size*size];
	protoMalloc(&in2D, size * size * sizeof(double));
	protoMalloc(&out2D, size * size * sizeof(double));

	protoLaunchKernel(initOne,1, NN, in2D,size*size);	
	

	dim3 grid((size+16-1)/16,(size+16-1)/16,1);
	dim3 block(16,16,1);

	protoLaunchKernel(lapl2D, grid, block, in2D, out2D, size);

	protoLaunchKernel(initOne,1, NN, in2D,size*size);	
	protoLaunchKernel(lapl2DShared, grid, block, in2D, out2D, size);

	protoLaunchKernel(initOne,1, NN, in2D,size*size);	
	protoLaunchKernel(lapl2DCoal, (size*size+256-1)/256, 256, in2D, out2D, size);
	
	protoMemcpy(host2D,out2D, (size*size)*sizeof(double), protoMemcpyDeviceToHost);

	good =true;
	for(int j = 1; j<size-1 ; j++)
	{
	for(int i = 1; i<size-1 ; i++)
		if(host2D[i+j*size] != 2)
		{
			good = false;
			std::cout << " error host2D["<<i<<"," << j<<"]:["<<i+j*size<<"] = " << host2D[i+j*size] << std::endl;
		        break;	
		}
		if(good == false) break;
	}
	if(good==true) std::cout << " stencil2D checked " << std::endl;
#endif
	return 0;
}
