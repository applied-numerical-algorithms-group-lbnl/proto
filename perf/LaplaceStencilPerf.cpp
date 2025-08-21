#include "../include/Proto.H"
#include "cuda.h"
#include "TestFunctions.H"
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

    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    int domainSize = layout.domain().box().size(0);
    double dx = 1.0/domainSize;
    LevelBoxData<double, 1, HOST> srcData(layout, Point::Ones());
    LevelBoxData<double, 1, HOST> dstData(layout, Point::Zeros());
    LevelBoxData<double, 1, HOST> slnData(layout, Point::Zeros());
    srcData.initialize(f_phi, dx, k, offset);
    slnData.initialize(f_Lphi, dx, k, offset);

    Stencil<double> L = Stencil<double>::Laplacian();

    double T = 0;
    for (int ii = 0; ii < numIter; ii++)
    {
        auto start = high_resolution_clock::now();
        for (auto iter : layout)
        {
            const auto& src = srcData[iter];
            auto& dst = dstData[iter];
            dst |= L(src);
        }
        auto end = high_resolution_clock::now();
        slnData.increment(dstData, -1.0);
        PROTO_ASSERT(slnData.absMax() < 1e-12, "Error: CPU invocation failed to get the right answer");
        auto elapsed = duration_cast<microseconds>(end - start);
        T += elapsed.count();
    }

    return T/numIter;
}

#ifdef PROTO_CUDA
double profileCuda(const DisjointBoxLayout& layout, int numIter = 10)
{

    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    int domainSize = layout.domain().box().size(0);
    double dx = 1.0/domainSize;
    LevelBoxData<double, 1, DEVICE> srcData(layout, Point::Ones());
    LevelBoxData<double, 1, DEVICE> dstData(layout, Point::Zeros());
    LevelBoxData<double, 1, DEVICE> slnData(layout, Point::Zeros());
    srcData.initialize(f_phi, dx, k, offset);
    slnData.initialize(f_Lphi, dx, k, offset);

    Stencil<double> L = Stencil<double>::Laplacian();

    double T = 0;
    for (int ii = 0; ii < numIter; ii++)
    {
        cudaEvent_t start, stop;

        CHECK_DEVICE_ERROR(cudaEventCreate(&start));
        CHECK_DEVICE_ERROR(cudaEventCreate(&stop));

        CHECK_DEVICE_ERROR(cudaEventRecord(start));
        for (auto iter : layout)
        {
            const auto& src = srcData[iter];
            auto& dst = dstData[iter];
            dst |= L(src);
        }
        CHECK_DEVICE_ERROR(cudaEventRecord(stop));
        CHECK_DEVICE_ERROR(cudaEventSynchronize(stop));
        CHECK_DEVICE_ERROR(cudaGetLastError());
        slnData.increment(dstData, -1.0);
        PROTO_ASSERT(slnData.absMax() < 1e-12, "Error: CPU invocation failed to get the right answer");
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

    int domainSize = 256;
    auto boxSize = 128*Point::Ones();
    int numIter = 1;

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
