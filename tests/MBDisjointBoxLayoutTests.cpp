#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"

using namespace Proto;
TEST(MBDisjointBoxLayout, Iteration) {
    int domainSize = 64;
    int boxSize = 16;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    
    std::set<std::pair<unsigned int, Box>> correctPatchSet;
    Box patches = Box::Cube(domainSize).coarsen(boxSizeVect);
    for (unsigned int bi = 0; bi < domain.size(); bi++)
    {
        for (auto p : patches)
        {
            Box b(p,p);
            b = b.refine(boxSizeVect);
            std::pair<unsigned int, Box> patch(bi, b);
            correctPatchSet.insert(patch);
        }
    }
    
    int N = pow(domainSize / boxSize, DIM)*domain.size();
    int n0 = N / numProc();
    int r =  N % numProc();
    if (procID() < r) { n0++; }
    std::set<std::pair<unsigned int, Box>> patchSet;
    for (auto index : layout)
    {
        Box b = layout[index];
        auto block = layout.block(index);
        std::pair<unsigned int, Box> patch(block, b);
        EXPECT_EQ(b.sizes(), boxSizeVect);
        EXPECT_EQ(patchSet.find(patch), patchSet.end());
        EXPECT_FALSE(correctPatchSet.find(patch) == correctPatchSet.end());
        for (auto pi : patchSet)
        {
            if (pi.first == block)
            {
                EXPECT_TRUE((pi.second & b).empty());
            }
        }
        patchSet.insert(patch);
    }
    EXPECT_EQ(patchSet.size(), n0);
}

TEST(MBDisjointBoxLayout, PatchConnectivity)
{
    int domainSize = 32;
    int boxSize = 16;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    
    Box domainBox = Box::Cube(domainSize / boxSize);
    Box adjX  = domainBox.adjacent(Point::X(), 1);
    Box adjY  = domainBox.adjacent(Point::Y(), 1);
    Box adjXY = domainBox.adjacent(Point::X() + Point::Y(), 1);
    Box edgX  = domainBox.edge(Point::X(), 1);
    Box edgY  = domainBox.edge(Point::Y(), 1);
    Box edgXY = domainBox.edge(Point::X() + Point::Y(), 1);
    auto CW = CoordPermutation::cw();
    auto CCW = CoordPermutation::ccw();
    auto R = CW*CW;
    for (auto i1 : layout)
    {
        for (auto i2 : layout)
        {
            auto b1 = layout.block(i1);
            auto b2 = layout.block(i2);
            auto p1 = layout.point(i1);
            auto p2 = layout.point(i2);
            Point q2 = Point::Zeros();

            Point p = layout.connectivity(i1,i2);
            if (b1 == b2)
            {
                q2 = p2;
            }
            else if (b2 == (b1 + 1) % layout.numBlocks())
            {
                q2 = CW.rotateCell(p2, edgY, adjX);
            }
            else if (b2 == (b1 + layout.numBlocks() - 1) % layout.numBlocks())
            {
                q2 = CCW.rotateCell(p2, edgX, adjY);
            }
            else {
                q2 = R.rotateCell(p2, edgXY, adjXY); 
            }

            if (Box::Kernel(1).shift(p1).contains(q2))
            {
                EXPECT_EQ(p, q2 - p1);
            } else {
                EXPECT_EQ(p, Point::Zeros());
            }
        }
    }
}

TEST(MBDisjointBoxLayout, Find) {
    int domainSize = 64;
    int boxSize = 16;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    
    for (auto iter : layout)
    {
        Point patch = layout.point(iter);
        auto block = layout.block(iter);
        auto index = layout.find(patch, block);
        EXPECT_EQ(index, iter);
    }
}

TEST(MBDisjointBoxLayout, Coarsen) {
    int domainSize = 64;
    int boxSize = 16;
    int refRatio = 2;
    int numBlocks = XPOINT_NUM_BLOCKS;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    
    std::vector<Point> refRatios(numBlocks, Point::Ones(refRatio));
    auto crseLayout = layout.coarsen(refRatios);

    auto crseDomain = buildXPoint(domainSize/refRatio);
    MBDisjointBoxLayout crseLayoutSoln(crseDomain, Point::Ones(boxSize/refRatio));

    EXPECT_EQ(crseLayout.numBoxes(), layout.numBoxes());
    EXPECT_TRUE(crseLayout.compatible(layout));
    EXPECT_TRUE(crseLayout.compatible(crseLayoutSoln));
    for (auto iter : crseLayout)
    {
        EXPECT_EQ(crseLayout[iter], crseLayoutSoln[iter]);
    }
}

TEST(MBDisjointBoxLayout, BlockBoundaries_XPoint)
{
    int domainSize = 8;
    int boxSize = 4;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Point patchDomainSize = Point::Ones(domainSize / boxSize);
    Box patchDomain(patchDomainSize);
    std::cout << "patchDomain: " << patchDomain << std::endl;

    for (auto iter_i : layout)
    {
        for (auto iter_j : layout)
        {
            auto pi = layout.point(iter_i);
            auto bi = layout.block(iter_i);
            auto pj = layout.point(iter_j);
            auto bj = layout.block(iter_j);
            
            if (bi == 0 && bi != bj)
            {
                auto dij = layout.connectivity(iter_i, iter_j);
                auto dji = layout.connectivity(iter_j, iter_i);
                if (dij == Point::Zeros())
                {
                    std::cout << "no connection between " << pi << " and " << pj << "(bj = " << bj << ")" << std::endl;
                    EXPECT_EQ(dji, Point::Zeros());
                } else {
                    std::cout << "\npi: " << pi << " ,bi: " << bi << " | pj: " << pj << ", " << bj << std::endl;
                    std::cout << "\tdij: " << dij << ", dji: " << dji << std::endl;
                }
            }

        }
    }
}
#if 0
#if DIM == 3
TEST(MBDisjointBoxLayout, BlockBoundaries_CubedSphere)
{
    int domainSize = 8;
    int boxSize = 4;
    int thickness = 8;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(boxSize);
    boxSizeVect[radialDir] = min(thickness, boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Point patchDomainSize = Point::Ones(domainSize / boxSize);
    patchDomainSize[radialDir] = thickness > boxSize ? thickness / boxSize : 1;
    Box patchDomain(patchDomainSize);
    std::cout << "patchDomain: " << patchDomain << std::endl;

    for (auto iter_i : layout)
    {
        for (auto iter_j : layout)
        {
            auto pi = layout.point(iter_i);
            auto bi = layout.block(iter_i);
            auto pj = layout.point(iter_j);
            auto bj = layout.block(iter_j);
            
            if (bi == 0 && bi != bj)
            {
                auto dij = layout.connectivity(iter_i, iter_j);
                auto dji = layout.connectivity(iter_j, iter_i);
                if (dij == Point::Zeros())
                {
                    std::cout << "no connection between " << pi << " and " << pj << "(bj = " << bj << ")" << std::endl;
                    EXPECT_EQ(dji, Point::Zeros());
                } else {
                    std::cout << "\npi: " << pi << " ,bi: " << bi << " | pj: " << pj << ", " << bj << std::endl;
                    std::cout << "\tdij: " << dij << ", dji: " << dji << std::endl;
                }
            }

        }
    }

}
#endif
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
