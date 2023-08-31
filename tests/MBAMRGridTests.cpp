#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;

#if 1
TEST(MBAMRGrid, Construction) {
    int domainSize = 64;
    int boxSize = 16;
    int numBlocks = 5;
    int numLevels = 3;
    auto domain = buildXPoint(domainSize, numBlocks);
    std::vector<Point> boxSizeVect(numBlocks, Point::Ones(boxSize));
    std::vector<Point> refRatios(numLevels-1, Point::Ones(2));

    MBAMRGrid grid(domain, boxSizeVect, refRatios);
    
    for (int li = 0; li < numLevels; li++)
    {
        auto& layout = grid.getLevel(li);
        std::cout << "====================================================" << std::endl;
        std::cout << "Level: " << li << std::endl;
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            auto box = layout[iter];
            std::cout << "\tblock: " << block << " | box: " << box << std::endl;
        }
    }
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
