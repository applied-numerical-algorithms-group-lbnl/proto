#include <gtest/gtest.h>
#include "Proto.H"
#include "TestFunctions.H"

using namespace Proto;

/*
template<typename T, unsigned int C>
bool compareBoxData(
        const BoxData<T, C, HOST>& a_src,
        const BoxData<T, C, HOST>& a_dst)
{
    auto B = a_src.box() & a_dst.box();
    for (auto pt : B)
    {
        for (int cc = 0; cc < C; cc++)
        {
            T diff = abs(a_src(pt, cc) - a_dst(pt, cc));
            if (diff > 1e-12)
            {
                return false;
            }
        }
    }
    return true;
}
*/
TEST(AMRData, Initialize) {
    int domainSize = 32;
    int numLevels = 3;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    double dx = 1.0/domainSize;
    Point refRatio = Point::Ones(2);
    Point boxSize = Point::Ones(16);
    auto grid = telescopingGrid(domainSize, numLevels, refRatio, boxSize);
    AMRData<double, 1, HOST> hostData(grid, Point::Ones());
    hostData.initialize(dx, f_phi, k, offset);
#ifdef PROTO_ACCEL
    AMRData<double, 1, DEVICE> deviData(grid, Point::Ones());
    deviData.initialize(dx, f_phi, k, offset);
#endif
    for (int lvl = 0; lvl < numLevels; lvl++)
    {
        double dx_lvl = dx / pow(refRatio[0], lvl);
        for (auto iter : grid[lvl])
        {
            auto& hostData_i = hostData[lvl][iter];
            int N = hostData_i.size();
            Box B = hostData_i.box();
            BoxData<double, 1, HOST> soln_i(B);
            forallInPlace_p(f_phi, soln_i, dx_lvl, k, offset);
            EXPECT_TRUE(compareBoxData(soln_i, hostData_i));
#ifdef PROTO_ACCEL
            BoxData<double, 1, HOST> tmpData_i(B);
            auto& deviData_i = deviData[lvl][iter];
            deviData_i.copyTo(tmpData_i);
            EXPECT_TRUE(compareBoxData(soln_i, tmpData_i));
#endif
        }
    }
}

TEST(AMRData, InitConvolve)
{
    int domainSize = 32;
    int numLevels = 3;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    int numIter = 3;
    Point refRatio = Point::Ones(2);
    Point boxSize = Point::Ones(16);
    double hostErr[numIter];
#ifdef PROTO_ACCEL
    double deviErr[numIter];
#endif
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = 1.0/domainSize;
        auto grid = telescopingGrid(domainSize, numLevels, refRatio, boxSize);
        AMRData<double, 1, HOST> hostData(grid, Point::Ones());
        AMRData<double, 1, HOST> soln(grid, Point::Ones());
        AMRData<double, 1, HOST> error(grid, Point::Ones());
        Operator::initConvolve(hostData, dx, f_phi, k, offset);
        soln.initialize(dx, f_phi_avg, k, offset);
        hostErr[nn] = 0;
#ifdef PROTO_ACCEL 
        AMRData<double, 1, DEVICE> deviData(grid, Point::Ones());
        Operator::initConvolve(deviData, dx, f_phi, k, offset);
        deviErr[nn] = 0;
#endif
        for (int lvl = 0; lvl < numLevels; lvl++)
        {
            for (auto iter : grid[lvl])
            {
                auto& hostData_i = hostData[lvl][iter];
                auto& soln_i = soln[lvl][iter];
                auto& error_i = error[lvl][iter];
                hostData_i.copyTo(error_i);
                error_i -= soln_i;
            }
            hostErr[nn] = std::max(hostErr[nn], error[lvl].absMax());
#ifdef PROTO_ACCEL
            for (auto iter : grid[lvl])
            {
                auto& deviData_i = deviData[lvl][iter];
                auto& soln_i = soln[lvl][iter];
                auto& error_i = error[lvl][iter];
                deviData_i.copyTo(error_i);
                error_i -= soln_i;
            }
            deviErr[nn] = std::max(deviErr[nn], error[lvl].absMax());
#endif
        }
        domainSize *= 2;
    }
    double rate = 4;
    for (int ii = 1; ii < numIter; ii++)
    {
        double hostRate_i = log(hostErr[ii-1]/hostErr[ii])/log(2.0);
        EXPECT_TRUE(abs(rate - hostRate_i) < 0.1);
#ifdef PROTO_ACCEL
        double deviRate_i = log(deviErr[ii-1]/deviErr[ii])/log(2.0);
        EXPECT_TRUE(abs(rate - deviRate_i) < 0.1);
#endif
    }
}

TEST(AMRData, Algebraic)
{
    int domainSize = 32;
    int numLevels = 3;
    Point offset(1,2,3,4,5,6);
    Point k(1,2,3,4,5,6);
    double dx = 1.0/domainSize;
    Point refRatio = Point::Ones(2);
    Point boxSize = Point::Ones(16);
    auto grid = telescopingGrid(domainSize, numLevels, refRatio, boxSize);
    AMRData<double, 1, HOST> hostDataAdd(grid, Point::Ones());
    AMRData<double, 1, HOST> hostDataMlt(grid, Point::Ones());
    hostDataAdd.initialize(dx, f_phi, k, offset);
    hostDataMlt.initialize(dx, f_phi, k, offset);
    hostDataAdd += 7;
    hostDataMlt *= 7;
#ifdef PROTO_ACCEL
    AMRData<double, 1, DEVICE> deviDataAdd(grid, Point::Ones());
    AMRData<double, 1, DEVICE> deviDataMlt(grid, Point::Ones());
    deviDataAdd.initialize(dx, f_phi, k, offset);
    deviDataMlt.initialize(dx, f_phi, k, offset);
    deviDataAdd += 7;
    deviDataMlt *= 7;
#endif
    for (int lvl = 0; lvl < numLevels; lvl++)
    {
        double dx_lvl = dx / pow(refRatio[0], lvl);
        for (auto iter : grid[lvl])
        {
            auto& hostDataAdd_i = hostDataAdd[lvl][iter];
            auto& hostDataMlt_i = hostDataMlt[lvl][iter];
            Box B = hostDataAdd_i.box();
            BoxData<double, 1, HOST> solnAdd_i(B);
            BoxData<double, 1, HOST> solnMlt_i(B);
            forallInPlace_p(f_phi, solnAdd_i, dx_lvl, k, offset);
            forallInPlace_p(f_phi, solnMlt_i, dx_lvl, k, offset);
            solnAdd_i += 7;
            solnMlt_i *= 7;
            EXPECT_TRUE(compareBoxData(solnAdd_i, hostDataAdd_i));
            EXPECT_TRUE(compareBoxData(solnMlt_i, hostDataMlt_i));
#ifdef PROTO_ACCEL
            BoxData<double, 1, HOST> tmpDataAdd_i(B);
            BoxData<double, 1, HOST> tmpDataMlt_i(B);
            auto& deviDataAdd_i = deviDataAdd[lvl][iter];
            auto& deviDataMlt_i = deviDataMlt[lvl][iter];
            deviDataAdd_i.copyTo(tmpDataAdd_i);
            deviDataMlt_i.copyTo(tmpDataMlt_i);
            EXPECT_TRUE(compareBoxData(solnAdd_i, tmpDataAdd_i));
            EXPECT_TRUE(compareBoxData(solnMlt_i, tmpDataMlt_i));
#endif
        }
    }
}

#define NUMCOMPS 1
PROTO_KERNEL_START
void f_advectionExactF(
        Point& a_pt,
        Var<double,NUMCOMPS>& a_U,
        double a_h,
        double a_time)
{
    double r0 = .125;
    for (int comp = 0; comp < NUMCOMPS; comp++)
    {
        double rsq = 0.;
        for (int dir = 0; dir < DIM ; dir++)
        {
            double xcen = fmod(.5 + a_time,1.); 
            double xdir = a_pt[dir]*a_h + .5*a_h;
            double del; 
            if (xcen > xdir) del = min(abs(xdir - xcen),abs(xdir - xcen + 1.));
            if (xcen <= xdir) del = min(abs(xdir - xcen),abs(xdir - xcen - 1.));
            rsq += pow(del,2);
        }
        double r = sqrt(rsq);
        if (r > r0)
        {
            a_U(comp) = 0.;
        }
        else
        {
            a_U(comp) = pow(cos(M_PI*r/r0/2),6);
        }
    }
}
PROTO_KERNEL_END(f_advectionExactF, f_advectionExact);

TEST(AMRData, CopyToRefinement)
{
constexpr int TIME_STEP = 0;
constexpr int PLOT_NUM = 0;
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    PR_TIMER_SETFILE("test_copyTo.time.table");
    double t0 = .74609375;
    int maxRefs = 1;
    Array<bool, DIM> periodicity;
    periodicity.fill(true);
    int domainSize = 16;
    Point boxSize = Point::Ones(4);
    Point refRatio = Point::Ones(2);
    double dx0 = 1.0/domainSize;
    Point ghost = Point::Ones(0);
    
    // construct fine grid LBD.
    Box domainBoxC = Box::Cube(domainSize);
    ProblemDomain domainC(domainBoxC, periodicity);
    DisjointBoxLayout layoutC(domainC, boxSize);
    
    LevelBoxData<double > dataC(layoutC, ghost);
    dataC.initialize(f_pointID);
    //Operator::initConvolve(dataC, f_advectionExact, dx0, t0);

    auto domainBoxF = domainBoxC.refine(refRatio);
    auto subdomainBoxF = Box::Cube(domainSize); //domainSize * 2 / 2
    ProblemDomain domainF(domainBoxF, periodicity);
    std::vector<Point> cfPatches;
    auto patchKernel = Box::Cube(8);
    h5.writeLevel(dx0,dataC,"DATA_0");
    int plotNum = 1;
    for (int jj = 0; jj < patchKernel.size(1); jj++)
    {
        for (int ii = 0; ii < patchKernel.size(0)/2; ii++)
        {
            int kk = jj*patchKernel.size(0)+ii;

            cfPatches.push_back(patchKernel[kk]);
            DisjointBoxLayout layoutF(domainF, cfPatches, boxSize);

            auto layoutCF = layoutF.coarsen(refRatio);

            LevelBoxData<double> dataCF(layoutCF, ghost);
            dataCF.setVal(1e6);
            dataC.copyTo(dataCF);
            pr_out() << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
            pr_out() << "NOW PRINTING " << kk+1 << "-PATCH DATASET" << std::endl;
            for (auto iter : dataCF)
            {
                dataCF[iter].printData();
            }
            h5.writeLevel(dx0,dataCF,"DATA_%i", plotNum);
            plotNum++;
        }
    }
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
