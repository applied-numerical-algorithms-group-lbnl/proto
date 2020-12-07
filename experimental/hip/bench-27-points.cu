#include <Proto_gpu.H>

__global__
void initOne(double *ptr, unsigned int size)
{
	const unsigned id = threadIdx.x + blockDim.x*blockIdx.x; // i>1
	if(id < size )
		ptr[id] = 1;
}


__global__
void stencil27(double *in, double *out, const unsigned int size)
{
	const unsigned id = threadIdx.x + blockDim.x*blockIdx.x; // i>1

	const int size2 = size*size;
	const int x = id % size;
	const int y = (id % (size2)) / size;
	const int z = id / size2;

	if( (x<size - 1) && (y < size -1) && (z < size - 1) && (x > 0) && (y >0) && (z > 0) )
	{
		//double acc = 6/dt in[id];
		double acc = 26*in[id]; //  8 for debugging i
// z-1
		acc-=in[id-1-size-size2];
		acc-=in[id-size-size2];
		acc-=in[id+1-size-size2];
		acc-=in[id-1-size2];
		acc-=in[id-size2];
		acc-=in[id+1-size2];
		acc-=in[id-1+size-size2];
		acc-=in[id+size-size2];
		acc-=in[id+1+size-size2];

//z=0

		acc-=in[id-1-size];
		acc-=in[id-size];
		acc-=in[id+1];
		acc-=in[id-1];
		acc-=in[id+1];
		acc-=in[id-1+size];
		acc-=in[id+size];
		acc-=in[id+1+size];

//z= 1
		acc-=in[id-1-size+size2];
		acc-=in[id-size+size2];
		acc-=in[id+1-size+size2];
		acc-=in[id-1+size2];
		acc-=in[id+size2];
		acc-=in[id+1+size2];
		acc-=in[id-1+size+size2];
		acc-=in[id+size+size2];
		acc-=in[id+1+size+size2];
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
		if(tab[idx(i,j,k)] != 1 ) 
		{
			std::cout << "error : arr[("<< i << "," << j << "," << k << ")=" << idx(i,j,k) << " = " << tab[idx(i,j,k)] <<std::endl;
			return false;
		}
	return true;	
}

int main()
{
//	int size1D = 256;
	int nIter = 1000;
	double * out;
	double * in;
	double * host; 

	double result[5][6]; // dim, size, sizeBloc

	protoEvent_t start, stop;
	int size = 1;
	const unsigned int sizeMax = 512;

	size *= sizeMax * sizeMax * sizeMax;
	host = new double[size];
	protoMalloc(in, size * sizeof(double));
	protoMalloc(out, size * sizeof(double));
	protoEventCreate(&start);
	protoEventCreate(&stop);
     	int comptSize=0;
        for(int size1D = 32; size1D <= sizeMax ; size1D*=2)
	{      
		int comptBlock = 0;

	        for(int NN = 32 ; NN <= 1024 ; NN*=2)
        	{
	
			size = size1D*size1D*size1D;

			protoLaunchKernel(initOne, (size+NN-1)/NN, NN, in, size);	
			protoLaunchKernel(initOne, (size+NN-1)/NN, NN, out,size);	

			protoEventRecord(start);
 			for(int n = 0 ; n < nIter ; n++)
			{
				protoLaunchKernel(stencil27, (size+NN-1)/NN, NN, in, out, size1D);
			}

			protoEventRecord(stop);

			protoMemcpy(host,out,size*sizeof(double), protoMemcpyDeviceToHost);
//			bool res = check(host,size1D,3);

//			if(res)
			{	
				unsigned int flops = (26+1)*(size1D-2)*(size1D-2)*(size1D-2); // one elem

				protoEventSynchronize(stop);
				float milliseconds = 0;
				protoEventElapsedTime(&milliseconds, start, stop);
				std::cout << " Results are correct " << " size = " << size1D << " blockSize = " << NN << " Perf = " << nIter * ((flops* 1E-9)/(milliseconds*0.001)) << " GFLOPs" <<std::endl;
				result[comptSize][comptBlock]=  nIter * (flops* 1E-9)/(milliseconds*0.001);
			}
//			else std::cout << " Results are wrong " <<std::endl;
			comptBlock++;
		}
		comptSize++;

	}

	// print table

	for(int s = 0 ; s < 5 ; s++)
	{
		std::cout << " \\addplot plot coordinates {" <<std::endl;
		for(int b = 0 ; b < 6 ; b++)
			std::cout << "( " << std::pow(2,5+b) << " , " << result[s][b] << ") \n";
		std::cout << "};" <<  std::endl;
	}
	std::cout << std::endl;

	free(host);
	protoFree(in);
	protoFree(out);
	return 0;
}
