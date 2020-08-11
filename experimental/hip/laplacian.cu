#include <Proto_gpu.H>

__global__
void initOne(double *ptr, unsigned int size)
{
	for(unsigned int i = threadIdx.x +blockDim.x*blockIdx.x ; i < size ; i+= gridDim.x*blockDim.x )
		ptr[i] = 1;
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

__global__
void lapl2D(double *in, double *out, const unsigned int size)
{
	const unsigned id = threadIdx.x + blockDim.x*blockIdx.x; // i>1

	const int x = id % size;
	const int y = id / size;

	if( x<size - 1 && y < size -1 && x > 0 && y >0)
	{
		double acc = 6*in[id];
		acc-= in[id-1]; 
		acc-= in[id+1]; 
		acc-= in[id + size]; 
		acc-= in[id - size]; 
		out[id] = acc; //4*in[i] - in[i-1] - in[i+1]; // 4 for checking
	}
}



__global__
void lapl3D(double *in, double *out, const unsigned int size)
{
	const unsigned id = threadIdx.x + blockDim.x*blockIdx.x; // i>1

	const int size2 = size*size;
	const int x = id % size;
	const int y = (id % (size2)) / size;
	const int z = id / size2;

	if( (x<size - 1) && (y < size -1) && (z < size - 1) && (x > 0) && (y >0) && (z > 0) )
	{
		double acc = 8*in[id];
		acc-= in[id-1]; 
		acc-= in[id+1]; 
		acc-= in[id + size]; 
		acc-= in[id + size2]; 
		acc-= in[id - size]; 
		acc-= in[id - size2]; 
		out[id] = acc; //4*in[i] - in[i-1] - in[i+1]; // 4 for checking
	}
}


const int NN=256;

bool check(double *tab, int size, int dim)
{
	auto idx = [dim,size](int x, int y, int z) -> int
	{
		int id = x;
		id += dim > 1 ? size * y : 0;
		id += dim > 2 ? size * size * z : 0;
		return id;
	};
	int dir[3];
	dir[0] = size;
	dir[1] = dim > 1 ? size : 3; 
	dir[2] = dim > 2 ? size : 3;
	for(int k = 1 ; k < dir[2] - 1 ; k++)
	for(int j = 1 ; j < dir[1] - 1 ; j++)
	for(int i = 1 ; i < dir[0] - 1 ; i++)
		if(tab[idx(i,j,k)] != 2 ) 
		{
			std::cout << "error : arr[("<< i << "," << j << "," << k << ")=" << idx(i,j,k) << " = " << tab[idx(i,j,k)] <<std::endl;
			return false;
		}
	return true;	
}

int main()
{
	const int dim=3;
	int size1D = 256;
	int size = 1;
	double * out[dim];
	double * in[dim];
	double * host[3]; 

	protoEvent_t start[dim], stop[dim];

	for(int i = 0 ; i < dim ; i++)
	{	
		size *= size1D;
		host[i] = new double[size];
		protoMalloc(&in[i], size * sizeof(double));
		protoMalloc(&out[i], size * sizeof(double));
		protoEventCreate(&start[i]);
		protoEventCreate(&stop[i]);
	}
	
	size = 1;

	for(int i = 0 ; i < dim ; i++)
	{	
		size *= size1D;
		protoLaunchKernel(initOne,1, NN, in[i],size);	

		protoEventRecord(start[i]);
		if(i==0)
		protoLaunchKernel(lapl1DBasicV2, (size+NN-1)/NN, NN, in[i], out[i], size1D);

		if(i==1)
		protoLaunchKernel(lapl2D, (size+NN-1)/NN, NN, in[i], out[i], size1D);

		if(i==2)
		protoLaunchKernel(lapl3D, (size+NN-1)/NN, NN, in[i], out[i], size1D);
		protoEventRecord(stop[i]);

		protoMemcpy(host[i],out[i],size*sizeof(double), protoMemcpyDeviceToHost);
		

		bool res = check(host[i],size1D,i+1);

		if(res)
		{	
			int flops = 2*(i+1)+1; // one elem
			for(int j = 0 ; j < i+1; j++)
				flops *= size1D-2; // all elems

			protoEventSynchronize(stop[i]);
			float milliseconds = 0;
			protoEventElapsedTime(&milliseconds, start[i], stop[i]);
			std::cout << " Results are correct for Dim = " << i+1 << " Perf = " << (flops* 1E-9)/(milliseconds*0.001) << " GFLOPs" <<std::endl;
		}
		else std::cout << " Results are wrong for Dim = " << i+1 <<std::endl;
	}

	for(int i = 0 ; i < dim ; i++)
	{	
		free(host[i]);
		protoFree(in[i]);
		protoFree(out[i]);
	}

	return 0;
}
