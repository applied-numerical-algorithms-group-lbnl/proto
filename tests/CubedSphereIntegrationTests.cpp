#include <gtest/gtest.h>
#include "Proto.H"

using namespace Proto;

template<typename T, typename MAP, MemType MEM>
class ExponentialSource
{
public:
    ExponentialSource(const MBLevelMap<MAP, MEM>&, MBIndex) {}
    static constexpr unsigned int numState() { return 1; }
    static constexpr unsigned int numAux() { return 0; }
    static Point ghost() { return Point::Zeros(); }
    void init() {}
    void setTime(double) {}
    void operator()(BoxData<T, 1, MEM>& out, const BoxData<T, 1, MEM>& state, T scale)
    {
        state.copyTo(out);
        out *= scale;
    }
};

TEST(CubedSphereIntegration, GenericRK4UsesIntermediateStages)
{
    auto domain = CubedSphereShell::Domain(8, 8);
    MBDisjointBoxLayout layout(domain, Point::Ones(8));
    auto map = CubedSphereShell::Map(layout, Point::Zeros());
    const MBInterpOp interpolation;
    MBLevelRK4<ExponentialSource, MBMap_CubedSphereShell<HOST>, double> rk(map, interpolation);
    MBLevelBoxData<double, 1, HOST> state(layout, Point::Zeros());
    state.setVal(1.0);
    double dt = 0.1;
    rk.advance(state, dt);
    double expected = 1.0 + dt + dt*dt/2.0 + dt*dt*dt/6.0 + dt*dt*dt*dt/24.0;
    for (auto index : layout)
    {
        for (auto point : layout[index]) { EXPECT_NEAR(state[index](point), expected, 1e-14); }
    }
}

TEST(CubedSphereIntegration, ZeroWidthBlockBoundaries)
{
    auto domain = CubedSphereShell::Domain(12, 12);
    MBDisjointBoxLayout layout(domain, Point::Ones(12));
    Array<Point, DIM+1> ghost;
    ghost.fill(Point::Zeros());
    ghost[0] = Point::Ones(3);
    MBLevelBoxData<double, 1, HOST> data(layout, ghost);
    data.setVal(2.5);
    for (auto index : layout)
    {
        for (auto dir : Box::Kernel(1))
        {
            EXPECT_TRUE(data.bounds(index, dir).empty());
        }
    }
    data.exchange();
    for (auto index : layout)
    {
        for (auto point : layout[index]) { EXPECT_DOUBLE_EQ(data[index](point), 2.5); }
    }
}

TEST(CubedSphereIntegration, AngularInterpolationPreservesConstants)
{
    auto domain = CubedSphereShell::Domain(12, 12);
    MBDisjointBoxLayout layout(domain, Point::Ones(12));
    Point ghost = Point::Ones(4);
    auto interp = CubedSphereShell::MHDInterpOp<HOST>(layout, ghost);
    MBLevelBoxData<double, 1, HOST> source(layout, ghost);
    MBLevelBoxData<double, 1, HOST> dest(layout, ghost);
    source.setVal(2.5);
    source.exchange();
    dest.setVal(-1.0);
    interp.apply(dest, source);
    for (const auto& op : interp.getOps())
    {
        const auto& target = op->target();
        EXPECT_NEAR(dest[target.dstIndex()](target.dstPoint()), 2.5, 1e-11);
    }
}

TEST(CubedSphereIntegration, CoarsenAcrossDifferentPatchSizes)
{
    auto domain = CubedSphereShell::Domain(12, 12);
    MBDisjointBoxLayout fineLayout(domain, Point::Ones(6));
    auto coarseDomain = domain.coarsen(Point::Ones(2));
    MBDisjointBoxLayout coarseLayout(coarseDomain, Point::Ones(6));
    MBLevelBoxData<double, 1, HOST> fine(fineLayout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> coarse(coarseLayout, Point::Zeros());
    for (auto index : fineLayout)
    {
        for (auto point : fineLayout[index])
        {
            fine[index](point) = 100.0*fineLayout.block(index) + point.sum();
        }
    }
    fine.coarsenTo(coarse, Point::Ones(2));
    for (auto index : coarseLayout)
    {
        for (auto point : coarseLayout[index])
        {
            EXPECT_DOUBLE_EQ(coarse[index](point),
                100.0*coarseLayout.block(index) + 2.0*point.sum() + 0.5*DIM);
        }
    }
}

TEST(CubedSphereIntegration, CheckpointMetadataAndLegacyIteration)
{
    auto domain = CubedSphereShell::Domain(8, 8);
    MBDisjointBoxLayout layout(domain, Point::Ones(8));
    MBLevelBoxData<double, 1, HOST> original(layout, Point::Zeros());
    MBLevelBoxData<double, 1, HOST> restored(layout, Point::Zeros());
    original.setVal(3.25);
    HDF5Handler writer;
    writer.setTime(1.25);
    writer.setTimestep(0.125);
    writer.setIter(7);
    const std::string filename = "cubed_integration_checkpoint";
    writer.writeMBLevel({"state"}, original, filename);
    HDF5Handler reader;
    reader.readMBLevel(restored, filename);
    EXPECT_DOUBLE_EQ(reader.time(), 1.25);
    EXPECT_DOUBLE_EQ(reader.dt(), 0.125);
    EXPECT_EQ(reader.iter(), 7);
    for (auto index : layout)
    {
        for (auto point : layout[index]) { EXPECT_DOUBLE_EQ(restored[index](point), 3.25); }
    }
    barrier();
    if (procID() == 0)
    {
        auto file = H5Fopen((filename + ".hdf5").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        auto group = H5Gopen(file, "/level_0", H5P_DEFAULT);
        EXPECT_GE(H5Adelete(group, "iteration"), 0);
        H5Gclose(group);
        H5Fclose(file);
    }
    barrier();
    reader.readMBLevel(restored, filename);
    EXPECT_EQ(reader.iter(), 0);
    EXPECT_DOUBLE_EQ(reader.time(), 1.25);
}

TEST(CubedSphereIntegration, AnisotropicExchangePreservesBlockCoordinates)
{
    auto domain = CubedSphereShell::Domain(30, 60, 0);
    MBDisjointBoxLayout layout(domain, Point::Ones(15));
    for (int radialGhost : {0, 1, 7})
    {
        Point ghost = Point::Basis(0, radialGhost)
                    + Point::Basis(1, 5) + Point::Basis(2, 5);
        MBLevelBoxData<double, 1, HOST> data(layout, ghost);
        auto value = [](unsigned int block, Point point)
        {
            return 1000000.0*block + 10000.0*point[0] + 100.0*point[1] + point[2];
        };
        // Give physical ghosts known values too, so every exchanged cell can
        // be checked independently of physical boundary conditions.
        for (auto index : layout)
            for (auto point : data[index].box())
                data[index](point) = value(layout.block(index), point);
        data.exchange();
        for (auto index : layout)
            for (auto dir : Box::Kernel(1))
                for (const auto& bound : data.bounds(index, dir))
                {
                    auto adjBlock = layout.block(bound.adjIndex);
                    Point origin = bound.localData->box().low();
                    auto adjOrigin = domain.convertPoint(
                        origin, layout.block(index), adjBlock, PR_CELL);
                    Array<Point, DIM> axes;
                    for (int axis = 0; axis < DIM; ++axis)
                        axes[axis] = domain.convertPoint(origin + Point::Basis(axis),
                            layout.block(index), adjBlock, PR_CELL) - adjOrigin;
                    // Block coordinate transforms are affine. Cache the axes
                    // instead of repeating expensive graph queries per cell.
                    for (auto point : bound.localData->box())
                    {
                        Point adjPoint = adjOrigin;
                        for (int axis = 0; axis < DIM; ++axis)
                            adjPoint += axes[axis] * (point[axis] - origin[axis]);
                        EXPECT_DOUBLE_EQ((*bound.localData)(point), value(adjBlock, adjPoint));
                    }
                }
    }
}

int main(int argc, char* argv[])
{
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
