#include "../include/Proto.H"
#include "cuda.h"
#include <chrono>

using namespace Proto;

namespace {

    #ifdef PROTO_CUDA
    #define CHECK_DEVICE_ERROR(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUA error in %s:%d: %s\n", __FILE__, __LINE__, \
                cudaGetErrorString(err)); \
            exit(EXIT_FAILURE); \
        } \
    } while(0)
    #else
    #define CHECK_DEVICE_ERROR(call)
    #endif

    DisjointBoxLayout testLayout(int domainSize, Point boxSize, std::set<Point> skipPatches)
    {
        Box domainBox = Box::Cube(domainSize); 
        Box patchBox = domainBox.coarsen(boxSize);
        std::vector<Point> patches;
        for (auto patch : patchBox)
        {
            bool skip = skipPatches.find(patch) != skipPatches.end();
            if (!skip) { patches.push_back(patch); }
        }
        std::array<bool, DIM> periodicity;
        periodicity.fill(true);
        ProblemDomain domain(domainBox, periodicity);
        return DisjointBoxLayout(domain, patches, boxSize);
    }
}

double profileCPU(const DisjointBoxLayout& layout, int numIter = 10)
{
    using namespace chrono;

    LevelBoxData<double, 1, HOST> data(layout, Point::Zeros());
    for (auto iter : data)
    {
        data[iter].setVal(iter.global());
    }

    double T = 0;
    for (int ii = 0; ii < numIter; ii++)
    {
        auto start = high_resolution_clock::now();
        auto absMaxValue = data.absMax();
        auto end = high_resolution_clock::now();
        PROTO_ASSERT(absMaxValue == layout.numBoxes()-1, "Error: CPU invocation failed to get the right answer");
        auto elapsed = duration_cast<microseconds>(end - start);
        T += elapsed.count();
    }

    return T/numIter;
}

#ifdef PROTO_CUDA
double profileCuda(const DisjointBoxLayout& layout, int numIter = 10)
{
    LevelBoxData<double, 1, DEVICE> data(layout, Point::Zeros());
    for (auto iter : data)
    {
        data[iter].setVal(iter.global());
    }

    double T = 0;
    auto absMaxValue = data.absMax();  //prewarm
    for (int ii = 0; ii < numIter; ii++)
    {
        cudaEvent_t start, stop;

        CHECK_DEVICE_ERROR(cudaEventCreate(&start));
        CHECK_DEVICE_ERROR(cudaEventCreate(&stop));

        CHECK_DEVICE_ERROR(cudaEventRecord(start));

        absMaxValue = data.absMax();

        CHECK_DEVICE_ERROR(cudaEventRecord(stop));
        CHECK_DEVICE_ERROR(cudaEventSynchronize(stop));
        CHECK_DEVICE_ERROR(cudaGetLastError());
        PROTO_ASSERT(absMaxValue == layout.numBoxes()-1, "Error: CPU invocation failed to get the right answer");
        float milliseconds = 0;
        CHECK_DEVICE_ERROR(cudaEventElapsedTime(&milliseconds, start, stop));
        T += milliseconds;
    }


    return T*1e3/numIter; //microseconds
}
#endif

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif

    int domainSize = 512;
    auto boxSize = 256*Point::Ones();
    int numIter = 10;

    std::set<Point> skipPatches;
    auto layout = testLayout(domainSize, boxSize, skipPatches);

    auto cpuTime = profileCPU(layout, numIter);
    std::cout << "CPU Time: " << cpuTime << " us (averaged over " << numIter << " iterations)\n";
    #ifdef PROTO_CUDA
    auto cudaTime = profileCuda(layout);
    std::cout << "Cuda Time: " << cudaTime << " us (averaged over " << numIter << " iterations)\n";
    #endif

#ifdef PR_MPI
    MPI_Finalize();
#endif

}
