#include <cstdio>
#include <Proto_gpu.H>

//typedef int (*funcptr) ();

__device__ int f() { return 42 ; }

typedef decltype(&f) funcptr;

__device__ funcptr f_ptr = f ;

__global__ void kernel ( funcptr func )
{
    int k = func () ;
    printf ("%d\n", k) ;
}


int main ()
{
    funcptr h_funcptr ;

    if (protoSuccess != protoMemcpyFromSymbol (&h_funcptr, (const void *)f_ptr, sizeof (funcptr), 0, protoMemcpyDeviceToHost))
        printf ("FAILED to get SYMBOL\n");

    protoLaunchKernel(kernel, 1, 1, (h_funcptr)) ;
    if (protoDeviceSynchronize() != protoSuccess)
        printf ("FAILED\n");
    else
        printf ("SUCCEEDED\n");
}
