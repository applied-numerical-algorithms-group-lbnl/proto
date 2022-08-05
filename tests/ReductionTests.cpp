#include <gtest/gtest.h>
#include "Proto.H"
#include <iostream>
using namespace Proto;

constexpr unsigned int bufferSize = 8;
int* hostBuffer;
int sumValue;
int sumAbsValue;
int maxValue;
int minValue;
int absMaxValue;

#ifdef PROTO_CUDA
int* deviBuffer;
#endif
void initBuffers()
{
    int totalBufferSize = bufferSize*numProc();
    hostBuffer = (int*)proto_malloc<HOST>(bufferSize*sizeof(int));
    for (int ii = 0; ii < bufferSize; ii++)
    {
        hostBuffer[ii] = ii + procID()*bufferSize;
        if (ii % 2 == 0) { hostBuffer[ii] *= -1; }
    }
    sumValue = totalBufferSize / 2;
    sumAbsValue = totalBufferSize / 2 * (totalBufferSize - 1);
    maxValue = totalBufferSize - 1;
    minValue = -(maxValue - 1);
    absMaxValue = maxValue;

#ifdef PROTO_CUDA
    deviBuffer = (int*)proto_malloc<DEVICE>(bufferSize*sizeof(int));
    proto_memcpy<HOST, DEVICE>(hostBuffer, deviBuffer, bufferSize*sizeof(int));
#endif
}

template <MemType MEM>
void max_wrapper(int* buffer) {
    Reduction<int, Max, MEM> rxn;
    rxn.reduce(buffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, maxValue);
}

TEST(Reduction, DeviceMax) {
    max_wrapper<HOST>(hostBuffer);
#ifdef PROTO_CUDA
    max_wrapper<DEVICE>(deviBuffer);
#endif
}

template <MemType MEM>
void min_wrapper(int* buffer) {
    Reduction<int, Min, MEM> rxn;
    rxn.reduce(buffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, minValue);
}

TEST(Reduction, DeviceMin) {
    min_wrapper<HOST>(hostBuffer);
#ifdef PROTO_CUDA
    min_wrapper<DEVICE>(deviBuffer);
#endif
}

template <MemType MEM>
void abs_wrapper(int* buffer) {
    Reduction<int, Abs, MEM> rxn;
    rxn.reduce(buffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, absMaxValue);
}

TEST(Reduction, DeviceAbsMax) {
    abs_wrapper<HOST>(hostBuffer);
#ifdef PROTO_CUDA
    abs_wrapper<DEVICE>(deviBuffer);
#endif
}
template <MemType MEM>
void sum_wrapper(int* buffer) {
    Reduction<int, Sum, MEM> rxn;
    rxn.reduce(buffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, sumValue);
}

TEST(Reduction, DeviceSum) {
    sum_wrapper<HOST>(hostBuffer);
#ifdef PROTO_CUDA
    sum_wrapper<DEVICE>(deviBuffer);
#endif
}

template <MemType MEM>
void sumabs_wrapper(int* buffer) {
    Reduction<int, SumAbs, MEM> rxn;
    rxn.reduce(buffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, sumAbsValue);
}

TEST(Reduction, DeviceSumAbs) {
    sumabs_wrapper<HOST>(hostBuffer);
#ifdef PROTO_CUDA
    sumabs_wrapper<DEVICE>(deviBuffer);
#endif
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    initBuffers();
    int result = RUN_ALL_TESTS();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
