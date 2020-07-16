#include <Proto_gpu.H> //useless

#include<cstdio>
#include<iostream>
#include<string>

void printCPU(std::string message)
{
	std::cout << " CPU : " << message << std::endl;
}
__global__
void printGPU()
{
	printf("Hello World from GPU!\n");
	//printf(" GPU : %s \n", message);
	//std::cout << " GPU : " << message << std::endl;
}

int main()
{
	//cpu
	printCPU ( " Bonjour le monde ");
	//gpu
	//const char* ptr_h =" Bonjour le monde ";
	protoLaunchKernel(printGPU, 1, 1);
	return 0;
}
