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
void justAdd( int* value)
{
	*value += 1;
}

template<typename Func, typename... Args>
void Launch(Func& Ker, int nbBlocks, int nbThreads, Args...args)
{
	protoLaunchKernelMemAsync( Ker, nbBlocks, nbThreads, 0, 0, args...);
}

int main()
{
	//cpu
	printCPU ( " Test: justAdd ");
	//gpu
	//const char* ptr_h =" Bonjour le monde ";
	int val_h=6;
	int* val_d;

	protoMalloc(&val_d,1*sizeof(int));

	printCPU ( "Value before ");
	printCPU ( val_h );
	protoMemcpy(val_d, &val_h, 1*sizeof(int), protoMemcpyHostToDevice);
	protoLaunchKernel(justAdd, 1, 1, val_d);
	Launch(justAdd, 1, 1, val_d);
	protoMemcpy(&val_h, val_d, 1*sizeof(int), protoMemcpyDeviceToHost);
	printCPU ( "Value after ");
	printCPU ( val_h );
	return 0;
}
