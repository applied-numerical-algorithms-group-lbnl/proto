#include <gtest/gtest.h>
#include "Proto.H"

using namespace Proto;

namespace {

    class TestStruct {
        public:
        TestStruct(){}
        TestStruct(Point p, Box b, double v)
        : point(p), box(b), value(v) {}

        static TestStruct Get(int proc)
        {
            Point pi = Point::Ones(proc);
            Box bi = Box::Cube(proc);
            double vi = proc;
            return TestStruct(pi, bi, vi);
        }

        bool operator<(const TestStruct& rhs) const
        {
            if (value != rhs.value) { return value < rhs.value; }
            else { return point < rhs.point; }
            return false;
        }

        Point point;
        Box box;
        double value;
    };
}


TEST(DistributedSet, InsertAndExchangeBalanced) {
    
    DistributedSet<TestStruct> pointSet;
    int proc = procID();
    TestStruct local = TestStruct::Get(proc);
    pointSet.insert(local);
    pointSet.exchange();
    EXPECT_EQ(pointSet.data().size(), numProc());
    for (int pi = 0; pi < numProc(); pi++)
    {
        EXPECT_EQ(pointSet.data().count(TestStruct::Get(pi)), 1);
    }
}
TEST(DistributedSet, InsertAndExchangeUnbalanced) {
    
    DistributedSet<TestStruct> pointSet;
    int proc = procID();
    if (proc % 2 == 0)
    {
        TestStruct local = TestStruct::Get(proc);
        pointSet.insert(local);
    }
    pointSet.exchange();
    EXPECT_EQ(pointSet.data().size(), (numProc()+1)/2);
    for (int pi = 0; pi < numProc(); pi++)
    {
        if (pi % 2 == 0)
        {
            EXPECT_EQ(pointSet.data().count(TestStruct::Get(pi)), 1);
        } 
        
    }
}
TEST(DistributedSet, ClearLocalAndGlobal) {
    
    DistributedSet<TestStruct> pointSet;
    int proc = procID();
    TestStruct local = TestStruct::Get(proc);
    pointSet.insert(local);
    pointSet.exchange();
    
    EXPECT_EQ(pointSet.data().size(), numProc());
    pointSet.clearGlobal();
    EXPECT_EQ(pointSet.data().size(), 0);
    pointSet.exchange();
    EXPECT_EQ(pointSet.data().size(), numProc());
    pointSet.clearLocal();
    pointSet.clearGlobal();
    EXPECT_EQ(pointSet.data().size(), 0);
    pointSet.exchange();
    EXPECT_EQ(pointSet.data().size(), 0);

}

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
