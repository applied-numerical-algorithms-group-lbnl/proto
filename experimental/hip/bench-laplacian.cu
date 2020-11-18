#include <Proto_gpu.H>

__global__
void initOne(double *ptr, unsigned int size)
{
	const unsigned id = threadIdx.x + blockDim.x*blockIdx.x; // i>1
	if(id < size )
		ptr[id] = 1;
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
		//double acc = 6/dt in[id];
		double acc = 8*in[id]; //  8 for debugging 
		acc-= in[id-1]; 
		acc-= in[id+1]; 
		acc-= in[id + size]; 
		acc-= in[id + size2]; 
		acc-= in[id - size]; 
		acc-= in[id - size2]; 
		out[id] = acc; 
	}
}



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
//	int size1D = 256;
	int size = 1;
	double * out[dim];
	double * in[dim];
	double * host[3]; 

	double result[dim][5][6]; // dim, size, sizeBloc

	protoEvent_t start[dim], stop[dim];
	const unsigned int sizeMax = 512;

	for(int i = 0 ; i < dim ; i++)
	{	
		size *= sizeMax;
		host[i] = new double[size];
		protoMalloc(&in[i], size * sizeof(double));
		protoMalloc(&out[i], size * sizeof(double));
		protoEventCreate(&start[i]);
		protoEventCreate(&stop[i]);
	}
     	int comptSize=0;
        for(int size1D = 32; size1D <= sizeMax ; size1D*=2)
	{      
		int comptBlock = 0;

	        for(int NN = 32 ; NN <= 1024 ; NN*=2)
        	{
	
			size = 1;

			for(int i = 0 ; i < dim ; i++)
			{	
				size *= size1D;
				protoLaunchKernel(initOne, (size+NN-1)/NN, NN, in[i], size);	
				protoLaunchKernel(initOne, (size+NN-1)/NN, NN, out[i],size);	

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
					std::cout << " Results are correct for Dim = " << i+1 << " size = " << size1D << " blockSize = " << NN << " Perf = " << (flops* 1E-9)/(milliseconds*0.001) << " GFLOPs" <<std::endl;
					result[i][comptSize][comptBlock]=  (flops* 1E-9)/(milliseconds*0.001);
				}
				else std::cout << " Results are wrong for Dim = " << i+1 <<std::endl;
			}
			comptBlock++;
		}
		comptSize++;

	}

	// print table

  	for(int d = 0 ; d < dim ; d++)
	{
		std::cout << " Table Dim = " << d+1 << std::endl << std::endl;
		for(int s = 0 ; s < 5 ; s++)
		{
			for(int b = 0 ; b < 6 ; b++)
				std::cout << result[d][s][b] << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	for(int i = 0 ; i < dim ; i++)
	{	
		free(host[i]);
		protoFree(in[i]);
		protoFree(out[i]);
	}
	return 0;
}
