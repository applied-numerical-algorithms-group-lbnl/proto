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

double profileCPU(const DisjointBoxLayout& layout)
{
    using namespace chrono;

    LevelBoxData<double, 1, HOST> data(layout, Point::Zeros());
    for (auto iter : data)
    {
        data[iter].setVal(iter.local());
    }

    auto start = high_resolution_clock::now();
    auto absMaxValue = data.absMax();
    auto end = high_resolution_clock::now();
    auto elapsed = duration_cast<milliseconds>(end - start);

    return elapsed.count();
}

#ifdef PROTO_CUDA
double profileCuda(const DisjointBoxLayout& layout)
{
    LevelBoxData<double, 1, DEVICE> data(layout, Point::Zeros());
    for (auto iter : data)
    {
        data[iter].setVal(iter.local());
    }
    cudaEvent_t start, stop;

    CHECK_DEVICE_ERROR(cudaEventCreate(&start));
    CHECK_DEVICE_ERROR(cudaEventCreate(&stop));

    CHECK_DEVICE_ERROR(cudaEventRecord(start));

    auto absMaxValue = data.absMax();
    CHECK_DEVICE_ERROR(cudaGetLastError());

    CHECK_DEVICE_ERROR(cudaEventRecord(stop));
    CHECK_DEVICE_ERROR(cudaEventSynchronize(stop));

    double milliseconds = 0;
    CHECK_DEVICE_ERROR(cudaEventElapsedTime(&milliseconds, start, stop));

    return milliseconds;
}
#endif

int main(int argc, char** argv)
{
    int domainSize = 512;
    auto boxSize = 32*Point::Ones();
    std::set<Point> skipPatches;
    auto layout = testLayout(domainSize, boxSize, skipPatches);

    auto cpuTime = profileCPU(layout);
    std::cout << "CPU Time: " << std::endl;
    #ifdef PROTO_CUDA
    auto cudaTime = profileCuda(layout);
    std::cout << "Cuda Time: " << std::endl;
    #endif

}
