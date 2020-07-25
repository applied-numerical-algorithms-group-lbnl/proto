#include <iostream>
#include <mpi.h>
#include<Proto_gpu.H>
//#include <omp.h>


template<typename F, class... T>
void loopStream(protoStream_t *s, int n, int* shift, F f, T... arg)
{
        for(int i=0 ; i<n ; i++)
                f(s[i], shift[i], shift[i+1], arg...);
}


/*template<typename F, class... T>
void loopStreamOMP(protoStream_t *s, int n, int* shift, F f, T... arg)
{
#pragma omp parallel for schedule(dynamic)
        for(int i=0 ; i<n ; i++)
                f(s[i], shift[i], shift[i+1], arg...);
}*/


__global__
void uselessKernel(double* t1, int n)
{

        int tid = blockIdx.x * blockDim.x + threadIdx.x;
        if(tid < n)
        {
                for(int i = 0 ; i< 10000 ; i++)
                        t1[tid] += 1;
                for(int i = 0 ; i< 10000 ; i++)
                        t1[tid] -= 1;


        }
}

void send(protoStream_t *s, double *ptr, int *shift, int nMsg, int to, int tag)
{
        for(int itMsg=0 ; itMsg < nMsg ; itMsg++)
        {
                int size = shift[itMsg+1] - shift[itMsg];

                // check that the unpack is over
                protoStreamSynchronize(s[itMsg]);

                // send
                MPI_Send(ptr + shift[itMsg], size, MPI_DOUBLE, to, tag, MPI_COMM_WORLD);
        }
}

void recv(double *ptr, int *shift, int nMsg, int from, int tag, MPI_Request *req )
{
        for(int itMsg=0 ; itMsg < nMsg ; itMsg++)
        {
                int size = shift[itMsg+1] - shift[itMsg];


                // send
                MPI_Irecv(ptr + shift[itMsg], size, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &req[itMsg]);
        }
}



int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);


    if(argc != 2)
    {
        std::cout << " error nb Arg " << std::endl;
        return 0;
    }


    std::size_t size=  std::stoi(argv[1],0) / sizeof(double);

    if(size >100000000)
    {
          std::cout << " size is too big " << std::endl;
          return 0;
    }


    bool protoSet = true;
    if(protoSet)
    {
        protoSetDevice(world_rank);
        if(world_rank==0) std::cout << " Master rank: use protoSetDevice " << std::endl;
    }
    else if(world_rank==0) std::cout << " Master rank: disable protoSetDevice " << std::endl;

    int nDevice=666;
    int nbDevice=667;
    protoGetDevice(&nDevice);
    protoGetDeviceCount(&nbDevice);

    // proto stream stuff
    const int nbStream = 8;
    protoStream_t streams[nbStream];
    for (int i = 0; i < nbStream; i++) protoStreamCreate(&streams[i]);

    // compute shift for each stream
    int nbItemPerStream = size / nbStream;
    int shift[nbStream+1];
    shift[nbStream]=size;
    for (int i = 0; i < nbStream; i++) shift[i]=i*nbItemPerStream;



    std::cout << " rank = " << world_rank << ", num gpu = " << nDevice << ", there are " << nbDevice << " gpu(s) " << std::endl;



    std::size_t nbBytes = size*sizeof(double);

    if(world_rank==0)  std::cout << "Master rank: number of Bytes " << nbBytes << std::endl;


    double* h_ptr;
    double* h_buff = new double[size];
    double* h_recv = new double[size];

    double* d_ptr;
    double* d_buff;
    double* d_recv;

    protoHostAlloc((void**)&h_ptr,nbBytes);

    protoMalloc((void**)&d_ptr,nbBytes);
    protoMalloc((void**)&d_buff,nbBytes);
    protoMalloc((void**)&d_recv,nbBytes);

     double t1, t2;

    //fill
    for(int it=0;it<size;it++)
            h_ptr[it]=world_rank;


    // MPI stuff

    MPI_Request requestsend, requestrecv, reqrecv[nbStream];
    MPI_Status status, stat[nbStream];


    // old version
    //
    //

    const int nbExchange = 100;
    const int nbBits = sizeof(double);

    MPI_Barrier(MPI_COMM_WORLD);

    t1 = MPI_Wtime();

    for(int step=0; step < nbExchange; step++)
    {
        protoDeviceSynchronize();

        // copy to device
        loopStream(streams, nbStream, shift, [](protoStream_t s, int begin, int end, double* host, double* device, int nbBits) {
                                protoMemcpyAsync(device+begin,host+begin,(end-begin)*nbBits,protoMemcpyHostToDevice,s);
                        }, h_ptr, d_ptr, nbBits);

        // kernel

        loopStream(streams, nbStream, shift, [](protoStream_t s, int begin, int end, double* tab) {
                int nbThreads = 1024;
                double * shift_ptr=tab+begin;
                int n = end-begin;
                int blockSize = n/nbThreads;
                uselessKernel<<< blockSize , nbThreads, 0, s>>>(shift_ptr,n);
        }, d_ptr);

        // copy to host
        loopStream(streams, nbStream, shift, [](protoStream_t s, int begin, int end, double* host, double* device, int nbBits) {
                                protoMemcpyAsync(host+begin,device+begin,(end-begin)*nbBits,protoMemcpyDeviceToHost,s);
                        }, h_ptr, d_ptr, nbBits);


        // ++ pack ++
        loopStream(streams, nbStream, shift, [](protoStream_t s, int begin, int end, double *in, double *out) {
                                protoStreamSynchronize(s);
                                std::copy(in+begin, in+end, out+begin);
                                }, h_ptr, h_buff);


        // send and recieve
        int from= (world_rank-1+world_size) % world_size ;
        int to = (world_rank+1) % world_size ;
        int tag=666+step;

        MPI_Irecv(h_recv, size, MPI_DOUBLE, from, tag, MPI_COMM_WORLD, &requestrecv);
        MPI_Send(h_buff, size, MPI_DOUBLE,   to, tag, MPI_COMM_WORLD);

        MPI_Wait(&requestrecv,&status);

        // ++ unpack ++
       loopStream(streams, nbStream, shift, [](protoStream_t s, int begin, int end, double *in, double *out) {
                             std::copy(in+begin, in+end, out+begin);
                             }, h_recv, h_ptr);

    }

    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();


    // Check
    // h_recv should be equal to rank - nbExchange modulo world_size

    double val = (world_rank - nbExchange + world_size*nbExchange) % world_size;
    for(int i = 0 ; i < size ; i++)
            if(h_ptr[i] != val)
            {
                    std::cout << "problem rank: " << world_rank << " indice " << i << ": " << h_ptr[i] << " != " << val << std::endl;
                   break;
            }

    MPI_Barrier(MPI_COMM_WORLD);


    if(world_rank==0)
    {
            std::cout << "Master rank: Host To Host OK \n" ;
            std::cout << "Master rank: " << nbExchange << " * ( HtoD + DtoH + std::copy + MPI recv/send + std::copy ) " << t2-t1 << " second(s) \n";

    }

    protoDeviceSynchronize();

    // new version
    //

    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();

    loopStream(streams, nbStream, shift, [](protoStream_t s, int begin, int end, double* host, double* device, int nbBits) {
                                protoMemcpyAsync(device+begin,host+begin,(end-begin)*nbBits,protoMemcpyHostToDevice,s);
                        }, h_ptr, d_ptr, nbBits);


    for(int step=0; step < nbExchange; step++)
    {

        // uselessKernel
        loopStream(streams, nbStream, shift, [](protoStream_t s, int begin, int end, double* tab) {
                int nbThreads = 1024;
                double * shift_ptr=tab+begin;
                int n = end-begin;
                int blockSize = n/nbThreads;
                uselessKernel<<< blockSize , nbThreads, 0, s>>>(shift_ptr,n);
        }, d_ptr);

        int from= (world_rank-1 + world_size) % world_size ;
        int to = (world_rank+1) % world_size ;
        int tag=566+step;

        // ++ pack ++ async
        loopStream(streams, nbStream, shift, [](protoStream_t s, int begin, int end, double* out, double* in, int nbBits) {
                                protoMemcpyAsync(out+begin,in+begin,(end-begin)*nbBits,protoMemcpyDeviceToDevice, s);
                        }, d_buff, d_ptr, nbBits);

         // send and recv
        recv(         d_recv, shift, nbStream, from, tag, reqrecv);
        send(streams, d_buff, shift, nbStream,   to, tag);

        MPI_Waitall(nbStream, reqrecv, stat);

        // unpack ++ async
        loopStream(streams, nbStream, shift, [](protoStream_t s, int begin, int end, double* out, double* in, int nbBits) {
                                protoMemcpyAsync(out+begin,in+begin,(end-begin)*nbBits,protoMemcpyDeviceToDevice,s);
                        }, d_ptr, d_recv, nbBits);

    }

    loopStream(streams, nbStream, shift, [](protoStream_t s, int begin, int end, double* host, double* device, int nbBits) {
                                protoMemcpyAsync(host+begin,device+begin,(end-begin)*nbBits,protoMemcpyDeviceToHost,s);
                        }, h_ptr, d_ptr, nbBits);



    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();



    protoFree(d_recv);
    protoFree(d_buff);
    protoFree(d_ptr);

    protoDeviceSynchronize();

    val = (world_rank - 2*nbExchange + 2*world_size*nbExchange) % world_size;
    for(int i = 0 ; i < size ; i++)
            if(h_ptr[i] != val)
            {
                    std::cout << "problem rank: " << world_rank << " indice " << i << ": " << h_ptr[i] << " != " << val << std::endl;
                   break;
            }

    MPI_Barrier(MPI_COMM_WORLD);


    if(world_rank==0)
    {
            std::cout << "Master rank: Device To Device OK \n" ;
            std::cout << "Master rank: HtoD + " << nbExchange << " * ( DtoD + MPI recv/send + DtoD) + DtoH = " << t2-t1 << " second(s) \n";

    }

    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}

