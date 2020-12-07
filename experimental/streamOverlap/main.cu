const int N = 1 << 25;

#include <iostream>
#include <pthread.h>
#include <omp.h>
#include <thread>

#include <chrono>
#include <thread>
#include "Proto_gpu.H"

__global__ void kernel(double *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    for (int i = 0; i < 80; i ++) {
        x[tid] = sqrt(pow(3.14159,i));
    }
}

__global__ void kernel2(double *x, int n)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
	for (int j = tid; j < n; j+=blockDim.x * gridDim.x)
    for (int i = 0; i < 80; i ++) {
        x[j] = sqrt(pow(3.14159,i));
    }
}



void launch_kernel(protoStream_t *streams, double **d_data, int i, int N)
{
	protoLaunchKernelMemAsync(kernel, N/1024, 1024, 0, streams[i], d_data[i], N);
	return;
}


void launch_copyhtod(protoStream_t *streams, double **h_data,  double **d_data, int i, int N)
{
	protoMemcpyAsync(d_data[i], h_data[i], N, protoMemcpyHostToDevice,streams[i]);
	return;
}

void launch_copydtoh(protoStream_t *streams, double **h_data,  double **d_data, int i, int N)
{
	protoMemcpyAsync(h_data[i], d_data[i], N, protoMemcpyDeviceToHost,streams[i]);
	return;
}

int main()
{
    const int num_streams = 4;

    protoStream_t streams[num_streams];

    std::thread threads[num_streams];



    for (int i = 0; i < num_streams; i++) {
        protoStreamCreate(&streams[i]);
    }

    double** h_data=   new double*[num_streams];
    double* d_data[num_streams];

   for(int i=0 ; i< num_streams; i++) {

     h_data[i] = new double[N];
     for(int j=0; j<N; j++)
	h_data[i][j]=1.;

     protoMalloc(d_data[i], N * sizeof(double));
   }

    
    int nbLoop = 3;

std::cout << " Init " << std::endl;
	

	for(int f=0; f<nbLoop; f++)
	{
		#pragma omp parallel for
		for(int i=0 ; i< num_streams; i++)     
		{  
			protoMemcpyAsync(d_data[i], h_data[i], N, protoMemcpyHostToDevice,streams[i]);
		}

		#pragma omp parallel for
        	for(int i=0 ; i< num_streams; i++)  
		{
			protoLaunchKernelMemAsync(kernel, N/64, 64, 0, streams[i], d_data[i], N);
		}

		#pragma omp parallel for
 		for(int i=0 ; i< num_streams; i++)
		{
			protoMemcpyAsync(h_data[i], d_data[i], N, protoMemcpyDeviceToHost,streams[i]); 
		}
		
	}


	protoDeviceSynchronize();
	std::this_thread::sleep_for(std::chrono::milliseconds(200));

	for(int f=0; f<nbLoop; f++)
	{
		#pragma omp parallel for
        	for(int i=0 ; i< num_streams; i++)  
		{      
			protoMemcpyAsync(d_data[i], h_data[i], N, protoMemcpyHostToDevice,streams[i]);

			protoLaunchKernelMemAsync(kernel, N/64, 64, 0, streams[i], d_data[i], N);

			protoMemcpyAsync(h_data[i], d_data[i], N, protoMemcpyDeviceToHost,streams[i]); 
		}

		
	}

    protoDeviceReset();

    return 0;
}
