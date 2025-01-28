#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"

using namespace Proto;

TEST(MBLevelArgs, GetAndSet) {
    int domainSize = 16;
    int boxSize = 8;
    int numBlocks = 5;
    int ghostWidth = 2;

    auto map = XPointMap<HOST>(domainSize, boxSize, numBlocks, ghostWidth);
    auto layout = map->layout();

    MBLevelArgs<MBMap_XPointRigid, HOST> args(map);

    auto d0 = std::make_shared<MBLevelBoxData<double, 2, DEVICE>>(layout, Point::Ones());
    auto d1 = std::make_shared<MBLevelBoxData<int, 3, HOST>>(layout, Point::Ones(2));

    int c0 = 17;
    double c1 = 3.33;

    args.Set("d0", d0);
    args.Set("d1", d1);

    args.Set("c0", c0);
    args.Set("c1", c1);

    auto& d00 = args.Get<double, 2, DEVICE>("d0");
    auto& d11 = args.Get<int, 3, HOST>("d1");
    auto c00 = args.Get<int>("c0");
    auto c11 = args.Get<double>("c1");

    // EXPECT_EQ(d0, &d00);
    // EXPECT_EQ(d1, &d11);

    // EXPECT_EQ(c0, c00);
    // EXPECT_EQ(c1, c11);
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
