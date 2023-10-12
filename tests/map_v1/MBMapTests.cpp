#include <gtest/gtest.h>
#include "Proto.H"
#include "Lambdas.H"

using namespace Proto;

template<typename Func>
bool testMapBoundData(MBMap<Func>& a_map)
{
    auto& layout = a_map.map().layout();
    const auto& domain = layout.domain();
    auto& graph = domain.graph();

    auto ghost = a_map.map().ghost();
    for (auto iter : layout)
    {
        auto locBlock = layout.block(iter);
        auto locBoxInterior  = layout[iter];
#if PR_VERBOSE > 0
        auto& locX = a_map.map()[iter];
        auto& locJ = a_map.jacobian()[iter];
#endif
        for (auto dir : Box::Kernel(1))
        {
            auto ghostDir = ghost[dir.codim()];
            auto XBounds = a_map.map().bounds(iter, dir);
            for (auto XBound : XBounds)
            {
                auto adjBoxInterior = layout[XBound.adjIndex];
                auto adjBlock = layout.block(XBound.adjIndex);
                auto adjDir = graph.reverseDir(locBlock, adjBlock, dir);
                auto locAdjRot = graph.rotation(locBlock, dir, adjBlock);
                auto adjLocRot = graph.rotation(adjBlock, dir, locBlock);

                Box locBoundBox0 = locBoxInterior.adjacent(dir*ghostDir);
                Box adjBoundBox0 = adjBoxInterior.edge(adjDir*ghostDir);
                
                Box locBoundBox1 = XBound.localData->box();
                Box adjBoundBox1 = XBound.adjData->box();
                
                auto adjBoundX = a_map(adjBoundBox0, adjBlock, adjBlock);
                auto locBoundX = a_map(locBoundBox0, adjBlock, locBlock);

                auto& adjX = a_map.map()[XBound.adjIndex];

                //EXPECT_EQ(locBoundBox0, locBoundBox1);
                //EXPECT_EQ(adjBoundBox0, adjBoundBox1);
               
                auto adjBoundErr = adjBoundX - (*XBound.adjData);
                auto locBoundErr = locBoundX - (*XBound.localData);

                //EXPECT_LT(adjBoundErr.absMax(), 1e-10);
                //EXPECT_LT(locBoundErr.absMax(), 1e-10);

#if PR_VERBOSE > 0
                pout() << "DIR: " << dir << std::endl;
                pout() << "adjBlock: " << adjBlock << " | adjDir: " << adjDir << std::endl;
                pout() << "---------------------------------------------------------------" << std::endl;
                pout() << "local bound box: " << locBoundBox0 << std::endl;
                pout() << "adjacent bound box: " << adjBoundBox0 << std::endl;
                pout() << "---------------------------------------------------------------" << std::endl;
                pout() << "local bound data: " << locBoundX.box() << "(" << locBoundX.box().size() << " points)" << std::endl;
                pout() << "adjacent bound data: " << adjBoundX.box() << "(" << adjBoundX.box().size() << " points)" << std::endl;
                pout() << "---------------------------------------------------------------" << std::endl;
                pout() << "local bound box (init): " << locBoundBox1 << std::endl;
                pout() << "adjacent bound box (init): " << adjBoundBox1 << std::endl;
                pout() << "---------------------------------------------------------------" << std::endl;
                pout() << "Adjacent Data: " << std::endl;
                adjX.printData();
                pout() << "Bound Data (loc): " << std::endl;
                XBound.localData->printData();
                pout() << "Bound Data (adj): " << std::endl;
                XBound.adjData->printData();

#endif

            }
            auto JBounds = a_map.jacobian().bounds(iter, dir);
            for (auto JBound : JBounds)
            {
                auto adjBoxInterior = layout[JBound.adjIndex];
                auto adjBlock = layout.block(JBound.adjIndex);
                auto adjDir = graph.reverseDir(locBlock, adjBlock, dir);
                auto locAdjRot = graph.rotation(locBlock, dir, adjBlock);
                auto adjLocRot = graph.rotation(adjBlock, dir, locBlock);

                Box locBoundBox0 = locBoxInterior.adjacent(dir*ghostDir);
                Box adjBoundBox0 = adjBoxInterior.edge(adjDir*ghostDir);
                
                Box locBoundBox1 = JBound.localData->box();
                Box adjBoundBox1 = JBound.adjData->box();

                auto& adjJ = a_map.jacobian()[JBound.adjIndex];
#if PR_VERBOSE > 0
                pout() << "DIR: " << dir << std::endl;
                pout() << "adjBlock: " << adjBlock << " | adjDir: " << adjDir << std::endl;
                pout() << "---------------------------------------------------------------" << std::endl;
                pout() << "local bound box (computed): " << locBoundBox0 << std::endl;
                pout() << "adjacent bound box (computed): " << adjBoundBox0 << std::endl;
                pout() << "---------------------------------------------------------------" << std::endl;
                pout() << "local bound box (init): " << locBoundBox1 << std::endl;
                pout() << "adjacent bound box (init): " << adjBoundBox1 << std::endl;
                pout() << "---------------------------------------------------------------" << std::endl;
                pout() << "Adjacent Data: " << std::endl;
                adjJ.printData(4);
                pout() << "Bound Data (loc): " << std::endl;
                JBound.localData->printData(4);
                pout() << "Bound Data (adj): " << std::endl;
                JBound.adjData->printData(4);

#endif
            }
        }
    }
    return true;
}


#if DIM > 2
template<typename Func>
bool testCubeSphere(MBMap<Func>& a_map, Point a_domainSize, double a_r0, double a_r1)
{
    bool success = true;

    auto& map = a_map.map();
    auto& layout = map.layout();
    const auto& domain = layout.domain();
    const auto& graph = domain.graph();

    Box B0(a_domainSize + Point::Ones());
    Point x = Point::Basis(0);
    Point y = Point::Basis(1);
    Point z = Point::Basis(2);
    
    for (auto iter : layout)
    {
        auto block = layout.block(iter);
        auto& patch = map[iter];
        
        Box b_r0 = B0.edge(-z) & patch.box();
        Box b_r1 = B0.edge(z)  & patch.box();
        for (auto pi : b_r0)
        {
            double r0 = sqrt(patch(pi, 0)*patch(pi, 0)
                           + patch(pi, 1)*patch(pi, 1)
                           + patch(pi, 2)*patch(pi, 2));
            success &= (abs(r0 - a_r0) < 1e-12);
        }
        for (auto pi : b_r1)
        {
            double r1 = sqrt(patch(pi, 0)*patch(pi, 0)
                           + patch(pi, 1)*patch(pi, 1)
                           + patch(pi, 2)*patch(pi, 2));
            success &= (abs(r1 - a_r1) < 1e-12);
        }
    }

    for (int b0 = 2; b0 < 6; b0++)
    {
        auto& L0 = layout.blockLayout(b0);
        auto& D0 = map.blockData(b0);
        int b1 = b0+1;
        if (b1 == 6) {b1 = 2;}
        auto& L1 = layout.blockLayout(b1);
        auto& D1 = map.blockData(b1);
        for (auto i0 : L0)
        {
            for (auto i1 : L1)
            {
                Point p0 = L0.point(i0);
                Point p1 = L1.point(i1);
                if (!L0.patchDomain().box().edge(x).contains(p0)) {continue; }
                if (!L1.patchDomain().box().edge(-x).contains(p1)) {continue; }
                if (p0[1] != p1[1]) {continue; }
                if (p0[2] != p1[2]) {continue; }
                auto& d0 = D0[i0];
                auto& d1 = D1[i1];
                Box b = B0.edge(-x) & d1.box();
                BoxData<double, 3, HOST> tmp0(b), tmp1(b);
                d0.copyTo(tmp0, d0.box(), -x*a_domainSize);
                d1.copyTo(tmp1);
                tmp0 -= tmp1;
                success &= tmp0.absMax() < 1e-12;
            }
        }
    }

    for (auto iter : layout)
    {
        auto& patch = map[iter];
        auto block = layout.block(iter);
        for (auto pi : patch.box())
        {
            Array<double, 3> x = patch.array(pi);
            Array<double, 3> X;
            for (int xi = 0; xi < 3; xi++)
            {
                X[xi] = (1.0*pi[xi])/(1.0*a_domainSize[xi]);
            }
            auto x_map = a_map(X, block);
            auto x_err = x_map - x;
            double err = x_err.absMax();
            EXPECT_LT(err, 1e-12);
        }

        for (auto pi : patch.box().extrude(Point::Ones(), -1))
        {
            Array<double, 3> x = patch.array(pi + Point::Ones());
            auto x_map = a_map(pi, Point::Ones(), block);
            auto x_err = x_map - x;
            double err = x_err.absMax();
            EXPECT_LT(err, 1e-12);
        }
    }

    for (unsigned int bi = 0; bi < 6; bi++)
    {
        Point dir_ij = Point::Basis(0);
        unsigned int bj = graph.adjacent(bi, dir_ij);
        Box box_0 = domain.blockDomain(bi).box();
        Box box_i = box_0.shift(dir_ij*(a_domainSize / 2));
        Point dir_ji = graph.connectivity(bj, bi);
        Box box_j = box_0.shift(dir_ji*(a_domainSize / 2));

        auto data_i = a_map(box_i, bi, bi);
        auto data_j = a_map(box_i, bi, bj);

        auto soln_i = forall_p<double, 3>(CubedSphereMap, box_i, bi, domain.blockSize(bi), PR_NODE);
        BoxData<double, 3> soln_j(box_j);
        soln_i.copyTo(soln_j, graph.rotation(bi, dir_ij, bj));
       
        BoxData<double, 3> err_i(box_i);
        BoxData<double, 3> err_j(box_j);
        data_i.copyTo(err_i);
        err_i -= soln_i;
        data_j.copyTo(err_j);
        err_j -= soln_j;
       
        double ei = err_i.absMax();
        double ej = err_j.absMax();

        EXPECT_LT(ei, 1e-12);
        EXPECT_LT(ej, 1e-12);
    }
    HDF5Handler h5;
    for (unsigned int bi = 0; bi < 6; bi++)
    {
        Box B = domain.blockDomain(bi).box();
        auto data = a_map(B, bi, PR_CELL);
        double err = 0.0;
        for (auto pi : B)
        {
            Array<double, DIM> x = pi;
            x += 0.5;
            Array<double, DIM> dx = a_domainSize;
            x /= dx;
            auto soln = a_map(x, bi);
            for (int dir = 0; dir < DIM; dir++)
            {
                err = max(err, abs(data(pi, dir) - soln[dir]));
            }
        }
        EXPECT_LT(err, 1e-12);
    }

    auto& J = a_map.jacobian();
    for (auto iter : layout)
    {
        auto& Ji = J[iter];
        EXPECT_GT(Ji.min(), 0);
    }

    return success;
}
#endif
#if 0
TEST(MBMap, Identity) {
    HDF5Handler h5;
    int domainSize = 8;
    int boxSize = 4;
    auto domain = buildIdentity(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    Point boundGhost = Point::Ones();
   
    // requires C++17
    MBMap<IdentityMap_t> map(IdentityMap, layout, ghost, boundGhost);
    
    auto& J = map.jacobian();
    double dx = 1.0/domainSize;
    double dv = pow(dx,DIM);
    for (auto iter : layout)
    {
        auto& Ji = J[iter];
        // J should be 1/(dx^DIM) everywhere
        EXPECT_LT(abs(Ji.max() - Ji.min()), 1e-12);
        EXPECT_LT(Ji.max()*pow(dx,DIM) - 1, 1e-12);
        EXPECT_GT(Ji.min(), 0);
    }
    for (unsigned int bi = 0; bi < 1; bi++)
    {
        Box B = domain.blockDomain(bi).box();
        auto data = map(B, bi, PR_CELL);
        double err = 0.0;
        for (auto pi : B)
        {
            Array<double, DIM> x = pi;
            x += 0.5;
            x *= dx;
            auto soln = map(x, bi);
            for (int dir = 0; dir < DIM; dir++)
            {
                err = max(err, abs(data(pi, dir) - soln[dir]));
            }
        }
        EXPECT_LT(err, 1e-12);
    }

    h5.writeMBLevel({"x", "y", "z"}, map.map(), "IDENTITY.map");
    h5.writeMBLevel({"J"}, map.jacobian(), "IDENTITY");
}
#endif
#if 0
TEST(MBMap, XPoint) {
    HDF5Handler h5;
    int domainSize = 8;
    int boxSize = 4;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    Point boundGhost = Point::Ones();
   
    MBMap<XPointMapRigid_t> map(XPointMapRigid, layout, ghost, boundGhost);
    //auto map = buildMap(XPointMap, layout, ghost);
    auto& J = map.jacobian();
    double dx = 1.0/domainSize;
    double dv = pow(dx,DIM);
    for (auto iter : layout)
    {
        auto& Ji = J[iter];
        EXPECT_GT(Ji.min(), 0);
    }

    h5.writeMBLevel({"x", "y", "z"}, map.map(), "XPOINT.map");
    h5.writeMBLevel({"J"}, map.jacobian(), "XPOINT");
}
#endif
#if DIM > 2
#if 0
TEST(MBMap, CubeSphere) {
    constexpr int C = 1;
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 8;
    auto domain = buildCubeSphere(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    Point boundGhost = Point::Ones();
   
    MBMap<AltCubedSphereMap_t> map(AltCubedSphereMap, layout, ghost, boundGhost);

    h5.writeMBLevel({"x", "y", "z"}, map.map(), "CUBE_SPHERE.map");
    h5.writeMBLevel({"J"}, map.jacobian(), "CUBE_SPHERE");

    //EXPECT_TRUE(testCubeSphere(map, Point::Ones(domainSize), 1.0, 2.0));
}
#endif
#if 1

PROTO_KERNEL_START
void f_initRadius(Point& a_pt, Var<double, 1>& a_r)
{
    a_r(0) = 1.0 + a_pt[0]*0.2;
}

TEST(MBMap, ThinCubeSphere) {
    HDF5Handler h5;
    
    int domainSize = 8;
    int boxSize = 8;
    int thickness = 1;
    
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones() - Point::Z());
    Point boundGhost = Point::Ones();
    //Point boundGhost = Point::Zeros();
    
    int N = 1;
    auto AVG = Stencil<double>::AvgDown(Point(2,2,2));
    for (int nn = 0; nn < N; nn++)
    {
        double dv = 1.0/(domainSize*domainSize*thickness);
        auto domain = buildThinCubeSphere(domainSize, thickness);
        Point boxSizeVect = Point::Ones(boxSize);
        //boxSizeVect[2] = thickness;
        boxSizeVect[0] = thickness;
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        MBMap<ThinCubedSphereMap_t> map(ThinCubedSphereMap, layout, ghost, boundGhost);
        MBLevelBoxData<double, 1, HOST> JTest(layout, ghost);
        MBLevelBoxData<double, 1, HOST> RTest(layout, ghost);
        MBLevelBoxData<double, 2, HOST> testData(layout, ghost);
        
        for (auto iter : layout)
        {
            auto block = layout.block(iter);
            Box B = layout[iter].grow(3);
            auto& J_i = JTest[iter];
            auto& R_i = RTest[iter];
            Array<BoxData<double, DIM>, DIM> NT;
            Array<double, DIM> h = Point::Ones();
            for (int ii = 0; ii < DIM; ii++)
            {
                NT[ii].define(layout[iter].grow(3));
                h[ii] /= (double)boxSizeVect[ii];
            }
            BoxData<double, 1> radius = forall_p<double, 1>(f_initRadius, B);
            radius.copyTo(R_i);
            BoxData<double, DIM> coords(layout[iter].grow(PR_NODE));

            Operator::cubedSphereGeometry(coords, NT, J_i, radius, layout[iter], h[0], block);
        }


        testData.initialize(f_MBPointID); 
        //testMapBoundData(map);

#if PR_VERBOSE > 0
        h5.writeMBLevel({"x", "y", "z"}, map, map.map(), "THIN_CUBE_SPHERE_MAP");
        h5.writeMBLevel({"J"}, map, map.jacobian(), "THIN_CUBE_SPHERE_J");
        h5.writeMBLevel({"J"}, map, JTest, "THIN_CUBE_SPHERE_JTEST");
        h5.writeMBLevel({"R"}, map, RTest, "THIN_CUBE_SPHERE_RTEST");
        h5.writeMBLevel({"i", "j"}, map, testData, "THIN_CUBE_SPHERE_ID");
#endif
        domainSize *= 2;
        boxSize *= 2;
        thickness *= 2;
    }
}
#endif
#endif
#if 0
TEST(MBMap, InitializeWithMap)
{
    HDF5Handler h5;
    int domainSize = 128;
    int boxSize = 64;
#if DIM > 2
    domainSize = 8;
    boxSize = 4;
#endif
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    Point boundGhost = Point::Ones();
   
    MBMap<XPointMap_t> map(XPointMap, layout, ghost, boundGhost);

    Point k{1,2,3,4,5,6};
    Array<double, DIM> offset{1,1,1,1,1,1};
    MBLevelBoxData<double, 1, HOST, PR_CELL> hostData(layout, ghost);
    hostData.initialize(f_phiM, map, k, offset);

    h5.writeMBLevel({"x", "y", "z"}, map.map(), "MAP_INIT.map");
    h5.writeMBLevel({"phi"}, hostData, "MAP_INIT");
}
#endif
#if 0
TEST(MBMap, AnalyticOps) {
    HDF5Handler h5;
    int domainSize = 8;
    int boxSize = 4;
    auto domain = buildXPoint(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    Point boundGhost = Point::Ones();
   
    MBMap<XPointMapRigid_t> map(XPointMapRigid, layout, ghost, boundGhost);
    
    for (auto iter : layout)
    {
        const auto& x_patch = map.map()[iter];
        unsigned int block = layout.block(iter);
        for (auto pi : x_patch.box().grow(-1))
        {
            auto x_low  = map(pi, Point::Ones(-1), block);
            auto x_mid = map(pi, Point::Zeros(), block);
            Array<double, DIM> mid = pi;
            mid += 0.5;
            mid /= domainSize;
            auto x_high = map(pi, Point::Ones(), block);
            for (int dir = 0; dir < DIM; dir++)
            {
                EXPECT_LT(x_patch(pi, dir) - x_low[dir], 1e-12);   
                EXPECT_LT(abs(map(mid, block)[dir] - x_mid[dir]), 1e-12);
                EXPECT_LT((x_patch(pi+Point::Ones(), dir) - x_high[dir]), 1e-12);
                
            }
        }
    }
}
TEST(MBMap, AnalyticOps2) {
    int domainSize = 8;
    int boxSize = 8;
    HDF5Handler h5;

    auto domain = buildShear(domainSize);
    Point boxSizeVect = Point::Ones(boxSize);
    MBDisjointBoxLayout layout(domain, boxSizeVect);

    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Ones(4));
    ghost[0] = Point::Ones(1);
    Point boundGhost = Point::Ones();

    // initialize map
    MBMap<ShearMap_t> map(ShearMap, layout, ghost, boundGhost);

    Box B0 = Box::Cube(boxSize);
    Box B1 = B0.adjacent(Point::X()).grow(1);

   
    auto x00 = map(B0,0,0);
    auto x10 = map(B0,1,0);
    auto x01 = map(B0,0,1);
    auto x11 = map(B0,1,1);
    
    pout() << "B0: " << B0 << std::endl;
    pout() << "map(B0,0,0): " << std::endl;
    x00.printData();
    pout() << "map(B0,1,0): " << std::endl;
    x10.printData();
    pout() << "map(B0,0,1): " << std::endl;
    x01.printData();
    pout() << "map(B0,1,1): " << std::endl;
    x11.printData();

    auto y00 = map(B1,0,0);
    auto y10 = map(B1,1,0);
    auto y01 = map(B1,0,1);
    auto y11 = map(B1,1,1);
    
    pout() << "\nB1: " << B1 << std::endl;
    pout() << "map(B1,0,0): " << std::endl;
    y00.printData();
    pout() << "map(B1,1,0): " << std::endl;
    y10.printData();
    pout() << "map(B1,0,1): " << std::endl;
    y01.printData();
    pout() << "map(B1,1,1): " << std::endl;
    y11.printData();
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
