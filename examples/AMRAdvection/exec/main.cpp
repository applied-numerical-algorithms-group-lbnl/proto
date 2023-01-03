#include <iostream>
#include <iomanip>
#include "Proto.H"
#include "BoxOp_Advection.H"
#include "AMRRK4.H"
#include "InputParser.H"

using namespace Proto;

int TIME_STEP = 0;
int PLOT_NUM = 0;
PROTO_KERNEL_START
void
f_initializeF(
        const Point& a_pt,
        Var<double,NUMCOMPS>& a_U,
        const double& a_h,
        const double& a_time)
{
    for (int comp = 0; comp < NUMCOMPS; comp++)
    {
        a_U(comp) = 1.;
        for (int dir = 0; dir < DIM ; dir++)
        {
            a_U(comp)*=(-cos(2*M_PI*(a_pt[dir]*a_h - a_time) + 2*M_PI*a_h) 
                      + cos(2*M_PI*(a_pt[dir]*a_h - a_time)))/(2*M_PI*a_h);
        }
    }
}
PROTO_KERNEL_END(f_initializeF, f_initialize);

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
            a_U(comp) = 1.;
        }
        else
        {
            a_U(comp) = 1.+pow(cos(M_PI*r/r0/2),6);
        }
    }
}
PROTO_KERNEL_END(f_advectionExactF, f_advectionExact);

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    using Proto::pout; 
    typedef BoxOp_Advection<double> OP;

    int domainSize = 32;
    int boxSize = 32;
    int refRatio = 4;
    int maxTimesteps = 64;
    int numLevels = 3;
    int regridBufferSize = 2;
    double maxTime = 1.0;
    int outputInterval = 1;
    double t0 = 0.0;
    Array<bool, DIM> periodicity;
    periodicity.fill(true);

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("boxSize",         boxSize);
    args.add("refRatio",        refRatio);
    args.add("maxTimesteps",    maxTimesteps);
    args.add("maxTime",         maxTime);
    args.add("outputInterval",  outputInterval);
    args.add("periodic_x",      periodicity[0]);
    args.add("periodic_y",      periodicity[1]);
    args.add("numLevels",       numLevels);
    args.add("regridBufferSize",regridBufferSize);
    args.add("t0",              t0);
    args.parse(argc, argv);
    args.print();

    if (boxSize > domainSize)
    {
        pout() << "Input error: boxSize > domainSize. Forcing boxSize == domainSize." << std::endl;
        boxSize = domainSize;
    }
    PR_TIMER_SETFILE(to_string(domainSize) 
        + ".DIM=" + to_string(DIM) + ".numProc=" + to_string(numProc())
        + "AMRAdvection.time.table");

    Point boxSizeVect = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeVect);

    double physDomainSize = 1.0;
    Array<double, DIM> dx;
    dx.fill(physDomainSize / domainSize);
    double dt = 0.5 / domainSize;

    pout() << setw(50) << setfill('=') << "=" << std::endl;
    pout() << "Coarsest Level Parameters" << std::endl;
    pout() << setw(50) << setfill('-') << "-" << std::endl;
    pout() << "dx: " << dx[0] << " | dt: " << dt << std::endl;
    pout() << "domain: " << domainBox << std::endl;
    pout() << setw(50) << setfill('-') << "-" << std::endl << std::endl;

    std::vector<Point> refRatios(numLevels-1, Point::Ones(refRatio));
    AMRGrid grid(layout, refRatios, numLevels);

    double dxLevel = dx[0];
    Point bufferSize = Point::Ones(regridBufferSize);
    for (int lvl = 0; lvl < numLevels-1; lvl++)
    {
        LevelBoxData<double, NUMCOMPS> initData(grid[lvl], Point::Zeros());
        //initData.initConvolve(f_advectionExact, dxLevel, t0);
        Operator::initConvolve(initData, f_advectionExact, dxLevel, t0);
        LevelTagData tags(grid[lvl], bufferSize);
        for (auto iter = grid[lvl].begin(); iter.ok(); ++iter)
        {
            OP::generateTags(tags[*iter], initData[*iter]);
        }
        AMRGrid::buffer(tags, bufferSize);
        grid.regrid(tags, lvl);
        dxLevel /= refRatio;
    }
    
    for (int lvl = 2; lvl < numLevels; lvl++)
    {
        grid.enforceNesting2(lvl);
    }

    AMRData<double, NUMCOMPS> U(grid, OP::ghost());
    //U.initConvolve(dx[0], f_advectionExact, t0);
    Operator::initConvolve(U, dx[0], f_advectionExact, t0);
    U.averageDown();
    AMRRK4<OP, double, NUMCOMPS> advectionOp(U, dx);

    double sum0 = U[0].sum();
    pout() << "Initial level 0 conservation sum: " << sum0 << std::endl;

    double time = t0;
#ifdef PR_HDF5
    h5.setTime(time);
    h5.setTimestep(dt);
    h5.writeAMRData(dx, U, "U_N0");
#endif
    for (int k = 0; ((k < maxTimesteps) && (time < maxTime)); k++)
    {
        TIME_STEP = k;
        PLOT_NUM = 0;
        advectionOp.advance(dt);
        time += dt;
#ifdef PR_HDF5
        h5.setTime(time);
#endif
        if ((k+1) % outputInterval == 0)
        {
            if (Proto::procID() == 0)
            {
                std::cout << "n: " << k << " | time: " << time << std::endl;
            }
#ifdef PR_HDF5
            h5.writeAMRData(dx, U, "U_N%i", k+1);
#endif
        }
    }

    sum0 = U[0].sum();
    pout() << "Final level 0 conservation sum: " << sum0 << std::endl;

    AMRData<double, NUMCOMPS> USoln(U.grid(), Point::Zeros());
    //USoln.initConvolve(dx[0], f_advectionExact, time);
    Operator::initConvolve(USoln, dx[0], f_advectionExact, time);
#ifdef PR_HDF5
    h5.writeAMRData(dx, USoln, "USoln");
#endif

    LevelBoxData<double, NUMCOMPS> UError(U.grid()[0], Point::Zeros());
    U[0].copyTo(UError);
    UError.increment(USoln[0], -1);
#ifdef PR_HDF5
    h5.writeLevel(dx, UError, "UError");
#endif
    double error = UError.absMax();
    std::cout << "error: " << error << std::endl;
    error /= dt;
    std::cout << "error/dt: " << error << std::endl;

    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}

