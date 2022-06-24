#include "ProtoAMR.H"
#include "InputParser.H"
#include "BoxOp_Advection.H"

#define NUMCOMPS 1

using namespace Proto;

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
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    HDF5Handler h5;

    int domainSize = 32;
    double physDomainSize = 1.0;
    int boxSize = 32;
    int numIter = 3;
    int numLevels = 3;
    int refRatio = 4;
    std::array<bool, DIM> periodicity;
    periodicity.fill(true);
    int regridBufferSize = 2;
    double t1 = 0.05;

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("physDomainSize",  physDomainSize);
    args.add("boxSize",         boxSize);
    args.add("numLevels",       numLevels);
    args.add("numIter",         numIter);
    args.add("refRatio",        refRatio);
    args.add("periodic_x",      periodicity[0]);
    args.add("periodic_y",      periodicity[1]);
#if DIM > 2
    args.add("periodic_z",      periodicity[2]);
#endif
    args.add("regridBufferSize",regridBufferSize);
    args.add("t1", t1);
    args.parse(argc, argv);
    args.print();

    typedef BoxOp_Advection<double> OP;
    
    std::array<double, DIM> dx;
    dx.fill(physDomainSize / domainSize);
    double dt = 0.5 / domainSize;
    double t0 = 0.0;

    Point boxSizeVect = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeVect);

    std::vector<Point> refRatios(numLevels-1, Point::Ones(refRatio));
    AMRGrid grid(layout, refRatios, numLevels);

    double dxLevel = dx[0];
    Point bufferSize = Point::Ones(regridBufferSize);
    for (int lvl = 0; lvl < numLevels-1; lvl++)
    {
        LevelBoxData<double, NUMCOMPS> initData(grid[lvl], Point::Zeros());
        initData.initConvolve(f_advectionExact, dxLevel, t0);
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
    AMRData<double, NUMCOMPS> U0(grid, OP::ghost());
    U0.initConvolve(dx[0], f_advectionExact, t0);
    U0.averageDown();
    h5.writeAMRData(dx, U0, "U0");

    AMRData<double, NUMCOMPS> U1(grid, OP::ghost());
    U1.initConvolve(dx[0], f_advectionExact, t1);
    U1.averageDown();

    h5.writeAMRData(dx, U1, "U1");
    
    double dx_i = dx[0];
    for (int lvl = numLevels-2; lvl >=0; lvl--)
    {
        LevelTagData tags(grid[lvl], bufferSize);
        for (auto iter = grid[lvl].begin(); iter.ok(); ++iter)
        {
            OP::generateTags(tags[*iter], U1[lvl][*iter]);
        }
        h5.writeLevel(dx_i, tags, "TAGS_L%i_0", lvl);
        grid.addFinerTags(tags, lvl); 
        h5.writeLevel(dx_i, tags, "TAGS_L%i_1", lvl);
        grid.regrid(tags, lvl);
        dx_i /= refRatio;
    }

#ifdef PR_MPI
    MPI_Finalize();
#endif
    return 0;
}













