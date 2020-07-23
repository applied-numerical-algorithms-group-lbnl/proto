#include <Proto_gpu.H> //useless

#include<cstdio>
#include<iostream>
#include<string>

template<typename T>
void printCPU(T message)
{
	std::cout << " CPU : " << message << std::endl;
}
__global__
void justAdd( double * ptr, int N, const double value)
{
	for(int i = threadIdx.x + blockIdx.x * blockDim.x ; i < N ; i += blockDim.x * gridDim.x )
		ptr[i]+=value;
}

template<typename Func, typename... Args>
void Launch(Func& Ker, int nbBlocks, int nbThreads, Args...args)
{
	hipLaunchKernelGGL( Ker, nbBlocks, nbThreads, 0, 0, args...);
}

int main()
{
	printCPU ( "-----------------");
	printCPU ( "| Test: streams |");
	printCPU ( "-----------------");


	const int nbStream = 8;
	const int sizePerStream = 10000;

	const int size = nbStream * sizePerStream;
	const int nbBytes = size * sizeof(double*); 


	double *hostArray;
	double *deviceArray;

	// init streams
	protoStream_t Streams[nbStream];
	
	for(int itStream = 0 ; itStream <  nbStream ; itStream++)
		protoStreamCreate(&Streams[itStream]);	

	// init array
	hostArray = new double[size];	
	protoHostAlloc(&deviceArray, nbBytes);

	for(int i = 0 ; i < size ; i++)
		hostArray[i] = 5;

	//copy array to device
	for(int itStream = 0 ; itStream <  nbStream ; itStream++)
		protoMemcpyAsync(deviceArray, hostArray, nbBytes, protoMemcpyHostToDevice, Streams[itStream]);
	
	double *ptr = deviceArray;
	for(int itStream = 0 ; itStream <  nbStream ; itStream++)
	{
		const double a = itStream;
		protoLaunchKernelMemAsync(justAdd, 1, 256, 0, Streams[itStream], ptr, sizePerStream, a);
		ptr += sizePerStream;
	}

	for(int itStream = 0 ; itStream <  nbStream ; itStream++)
	protoMemcpyAsync(hostArray, deviceArray, nbBytes, protoMemcpyDeviceToHost, Streams[itStream]);

	bool success = true;

	for(int itStream = 0 ; itStream <  nbStream ; itStream++)
		for(int i = 0; i < sizePerStream; i++)
			if(hostArray[itStream*sizePerStream + i] != 5 + itStream)
				success = false;
	if(success)
		printCPU(" The results are correct" );
	else
		printCPU(" The results are wrong");

	//end
	protoFree(deviceArray);
	for(int itStream = 0 ; itStream <  nbStream ; itStream++)
                protoStreamDestroy(Streams[itStream]);

	return 0;
}

