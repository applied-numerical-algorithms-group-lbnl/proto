const int N = 1 << 25;

#include <iostream>
#include <pthread.h>
#include <omp.h>
#include <thread>

#include <chrono>
#include <thread>

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



void launch_kernel(cudaStream_t *streams, double **d_data, int i, int N)
{
	kernel<<<N/1024, 1024, 0, streams[i]>>>(d_data[i], N);
	//cudaStreamSynchronize(streams[i]);
	return;
}


void launch_copyhtod(cudaStream_t *streams, double **h_data,  double **d_data, int i, int N)
{
	cudaMemcpyAsync(d_data[i], h_data[i], N, cudaMemcpyHostToDevice,streams[i]);
	//cudaStreamSynchronize(streams[i]);
	return;
}

void launch_copydtoh(cudaStream_t *streams, double **h_data,  double **d_data, int i, int N)
{
	cudaMemcpyAsync(h_data[i], d_data[i], N, cudaMemcpyDeviceToHost,streams[i]);
	//cudaStreamSynchronize(streams[i]);
	return;
}

int main()
{
    const int num_streams = 4;

    cudaStream_t streams[num_streams];

    std::thread threads[num_streams];



    for (int i = 0; i < num_streams; i++) {
        cudaStreamCreate(&streams[i]);
    }

    double** h_data=   new double*[num_streams];
    double* d_data[num_streams];

   for(int i=0 ; i< num_streams; i++) {

     h_data[i] = new double[N];
     for(int j=0; j<N; j++)
	h_data[i][j]=1.;

     cudaMalloc(&d_data[i], N * sizeof(double));
   }

    
    int nbLoop = 3;

std::cout << " Init " << std::endl;
	

	for(int f=0; f<nbLoop; f++)
	{
		/*for(int i=0 ; i< num_streams; i++)  
		{     
			threads[i] = std::thread(&launch_copyhtod, streams, h_data, d_data, i, N);
			threads[i].join();
		}


        	for(int i=0 ; i< num_streams; i++)  
		{
			//kernel<<<N/1024, 1024, 0, streams[i]>>>(d_data[i], N);
			threads[i] = std::thread(&launch_kernel, streams, d_data, i, N);
			threads[i].join();
		}

		for(int i=0 ; i< num_streams; i++)  
		{     
			threads[i] = std::thread(&launch_copydtoh, streams, h_data, d_data, i, N);
			threads[i].join();
		}*/

		#pragma omp parallel for
		for(int i=0 ; i< num_streams; i++)     
		{  
			
			cudaMemcpyAsync(d_data[i], h_data[i], N, cudaMemcpyHostToDevice,streams[i]);
			//cudaStreamSynchronize(streams[i]);
			
		}

		#pragma omp parallel for
        	for(int i=0 ; i< num_streams; i++)  
		{
		//	int ff=1024/num_streams;
		//	kernel2<<<1, ff, 0, streams[i]>>>(d_data[i], N);
			kernel<<<N/64, 64, 0, streams[i]>>>(d_data[i], N);
			//cudaStreamSynchronize(streams[i]);
			
		}

		#pragma omp parallel for
 		for(int i=0 ; i< num_streams; i++)
		{
			//cudaStreamSynchronize(streams[i]);
			cudaMemcpyAsync(h_data[i], d_data[i], N, cudaMemcpyDeviceToHost,streams[i]); 
			//cudaStreamSynchronize(streams[i]);
		}


		
	}


	cudaDeviceSynchronize();
std::this_thread::sleep_for(std::chrono::milliseconds(200));







	for(int f=0; f<nbLoop; f++)
	{
		#pragma omp parallel for
        	for(int i=0 ; i< num_streams; i++)  
		{      
			cudaMemcpyAsync(d_data[i], h_data[i], N, cudaMemcpyHostToDevice,streams[i]);

			kernel<<<N/64, 64, 0, streams[i]>>>(d_data[i], N);

			cudaMemcpyAsync(h_data[i], d_data[i], N, cudaMemcpyDeviceToHost,streams[i]); 
		}

		
	}

    cudaDeviceReset();

    return 0;
}
