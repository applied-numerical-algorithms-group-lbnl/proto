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

TEST(Reduction, HostMax) {
    initBuffers();
    Reduction<int, Max, HOST> rxn;
    rxn.reduce(hostBuffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, maxValue);
}

TEST(Reduction, HostMin) {
    Reduction<int, Min, HOST> rxn;
    rxn.reduce(hostBuffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, minValue);
}

TEST(Reduction, HostAbsMax) {
    Reduction<int, Abs, HOST> rxn;
    rxn.reduce(hostBuffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, absMaxValue);
}

TEST(Reduction, HostSum) {
    Reduction<int, Sum, HOST> rxn;
    rxn.reduce(hostBuffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, sumValue);
}

TEST(Reduction, HostSumAbs) {
    Reduction<int, SumAbs, HOST> rxn;
    rxn.reduce(hostBuffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, sumAbsValue);
}

#ifdef PROTO_CUDA

TEST(Reduction, DeviceMax) {
    initBuffers();
    Reduction<int, Max, DEVICE> rxn;
    rxn.reduce(deviBuffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, maxValue);
}

TEST(Reduction, DeviceMin) {
    Reduction<int, Min, DEVICE> rxn;
    rxn.reduce(deviBuffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, minValue);
}

TEST(Reduction, DeviceAbsMax) {
    Reduction<int, Abs, DEVICE> rxn;
    rxn.reduce(deviBuffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, absMaxValue);
}

TEST(Reduction, DeviceSum) {
    Reduction<int, Sum, DEVICE> rxn;
    rxn.reduce(deviBuffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, sumValue);
}

TEST(Reduction, DeviceSumAbs) {
    Reduction<int, SumAbs, DEVICE> rxn;
    rxn.reduce(deviBuffer, bufferSize);
    int result = rxn.fetch();
    EXPECT_EQ(result, sumAbsValue);
}

#endif

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    int result = RUN_ALL_TESTS();
#ifdef PR_MPI
    MPI_Finalize();
#endif
    return result;
}
