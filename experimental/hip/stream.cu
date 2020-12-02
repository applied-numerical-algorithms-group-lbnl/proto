/*
  STREAM benchmark implementation in CUDA.

    COPY:       a(i) = b(i)                 
    SCALE:      a(i) = q*b(i)               
    SUM:        a(i) = b(i) + c(i)          
    TRIAD:      a(i) = b(i) + q*c(i)        

  It measures the memory system on the device.
  The implementation is in single precision.

  Code based on the code developed by John D. McCalpin
  http://www.cs.virginia.edu/stream/FTP/Code/stream.c

  Written by: Massimiliano Fatica, NVIDIA Corporation
*/

#include <Proto_gpu.H> 
#define N	2000000
#define NTIMES	10


#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <sys/time.h>

# ifndef MIN
# define MIN(x,y) ((x)<(y)?(x):(y))
# endif
# ifndef MAX
# define MAX(x,y) ((x)>(y)?(x):(y))
# endif

static double	avgtime[4] = {0}, maxtime[4] = {0},
		mintime[4] = {FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};


static char	*label[4] = {"Copy:      ", "Scale:     ", "Add:       ", "Triad:     "};

static double	bytes[4] = {
    2 * sizeof(float) * N,
    2 * sizeof(float) * N,
    3 * sizeof(float) * N,
    3 * sizeof(float) * N
    };

/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.  */


double mysecond()
{
        struct timeval tp;
        struct timezone tzp;
        int i = gettimeofday(&tp,&tzp);
        return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


__global__ void set_array(float *a,  float value, int len)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < len) a[idx] = value;
}

__global__ void STREAM_Copy(float *a, float *b, int len)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < len) b[idx] = a[idx];
}

__global__ void STREAM_Scale(float *a, float *b, float scale,  int len)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < len) b[idx] = scale* a[idx];
}

__global__ void STREAM_Add( float *a, float *b, float *c,  int len)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < len) c[idx] = a[idx]+b[idx];
}

__global__ void STREAM_Triad( float *a, float *b, float *c, float scalar, int len)
{
  int idx = threadIdx.x + blockIdx.x * blockDim.x;
  if (idx < len) c[idx] = a[idx]+scalar*b[idx];
}

int main()
{
 float *d_a, *d_b, *d_c;
 int j,k;
 double times[4][NTIMES];
 float scalar;

 printf(" STREAM Benchmark implementation in CUDA\n");
 printf(" Array size (single precision)=%d\n",N);
 
 /* Allocate memory on device */
 protoMalloc((void**)&d_a, sizeof(float)*N);
 protoMalloc((void**)&d_b, sizeof(float)*N);
 protoMalloc((void**)&d_c, sizeof(float)*N);

 /* Compute execution configuration */
  dim3 dimBlock(128);
  dim3 dimGrid(N/dimBlock.x );
  if( N % dimBlock.x != 0 ) dimGrid.x+=1;

 printf(" using %d threads per block, %d blocks\n",dimBlock.x,dimGrid.x);

 /* Initialize memory on the device */
  protoLaunchKernel(set_array,dimGrid,dimBlock,d_a, 2.f, N);
  protoLaunchKernel(set_array,dimGrid,dimBlock,d_b, .5f, N);
  protoLaunchKernel(set_array,dimGrid,dimBlock,d_c, .5f, N);

/*	--- MAIN LOOP --- repeat test cases NTIMES times --- */

 scalar=3.0f;
 for (k=0; k<NTIMES; k++)
 {
  times[0][k]= mysecond();
  protoLaunchKernel(STREAM_Copy,dimGrid,dimBlock,d_a, d_c, N);
  protoThreadSynchronize();
  times[0][k]= mysecond() -  times[0][k];

  times[1][k]= mysecond();
  protoLaunchKernel(STREAM_Scale,dimGrid,dimBlock,d_b, d_c, scalar,  N);
  protoThreadSynchronize();
  times[1][k]= mysecond() -  times[1][k];

  times[2][k]= mysecond();
  protoLaunchKernel(STREAM_Add,dimGrid,dimBlock,d_a, d_b, d_c,  N);
  protoThreadSynchronize();
  times[2][k]= mysecond() -  times[2][k];

  times[3][k]= mysecond();
  protoLaunchKernel(STREAM_Triad,dimGrid,dimBlock,d_b, d_c, d_a, scalar,  N);
  protoThreadSynchronize();
  times[3][k]= mysecond() -  times[3][k];
 }

/*	--- SUMMARY --- */

    for (k=1; k<NTIMES; k++) /* note -- skip first iteration */
	{
	for (j=0; j<4; j++)
	    {
	    avgtime[j] = avgtime[j] + times[j][k];
	    mintime[j] = MIN(mintime[j], times[j][k]);
	    maxtime[j] = MAX(maxtime[j], times[j][k]);
	    }
	}

printf("Function      Rate (MB/s)   Avg time     Min time     Max time\n");
    for (j=0; j<4; j++) {
	avgtime[j] = avgtime[j]/(double)(NTIMES-1);

	printf("%s%11.4f  %11.4f  %11.4f  %11.4f\n", label[j],
	       1.0E-06 * bytes[j]/mintime[j],
	       avgtime[j],
	       mintime[j],
	       maxtime[j]);
    }


 /* Free memory on device */
 protoFree(d_a);
 protoFree(d_b);
 protoFree(d_c);

}
