#include "Proto.H"

using namespace Proto;

void f_foo(Point& a_pt, Var<double>& a_data, double a_dx, double a_k)
{
    double x[DIM];
    double r = 0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = a_pt[dir]*a_dx + a_dx / 2.0;
        r += x[dir];
    }
    a_data(0) = sin(2.0*M_PI*a_k*r);
}

int main(int argc, char** argv)
{
    // GENERAL PARAMETERS
    int domainSize = 64;
    int boxSize = 8;
    int refRatio = 2;
    double physDomainSize = 1.0;
    double k = 1.0;
    double dx[DIM];
    dx[0] = physDomainSize / domainSize;
    dx[1] = dx[0] / refRatio;

    // INTERP STENCIL PARAMETERS
    int interpStencilMaxShift = 2;
    Box interpStencilShiftKernel = Box::Kernel(1);
    int interpStencilPolynomialOrder = 2;

    // BUILD INTERPOLATION OPERATOR
    auto interpolate = InterpStencil<double>::Build(
        interpStencilMaxShift,
        interpStencilShiftKernel,
        interpStencilPolynomialOrder,
        refRatio);

    // BUILD GEOMETRY
    Point boxSizeVect = Point::Ones(boxSize);
    Point refRatioVect = Point::Ones(refRatio);
    Box domainBox = Box::Cube(domainSize);
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }
    Box refinedRegion = Box::Cube(domainSize / 2).shift(Point::Ones(domainSize / 4)).refine(refRatio);
    Box refinedPatches = refinedRegion.coarsen(boxSize);
    std::vector<Point> refinedPatchPoints;
    for (auto iter = refinedPatches.begin(); *iter != iter.end(); ++iter)
    {
        refinedPatchPoints.push_back(*iter);
    }

    ProblemDomain crseDomain(domainBox, periodicity);
    ProblemDomain fineDomain = crseDomain.refine(refRatio);
    DisjointBoxLayout crseLayout(crseDomain, boxSizeVect);
    DisjointBoxLayout fineLayout(fineDomain, refinedPatchPoints, boxSizeVect);
    DisjointBoxLayout crseFineLayout = fineLayout.coarsen(refRatioVect);    
    
    // INITIALIZE DATA HOLDERS
    LevelBoxData<double> crseData(crseLayout, Point::Ones());
    LevelBoxData<double> tempData(crseFineLayout, Point::Ones() + interpolate.spanPoint());
    LevelBoxData<double> fineData(fineLayout, Point::Ones());
    LevelBoxData<double> fineSoln(fineLayout, Point::Ones());
    
    crseData.initialize(f_foo, dx[0], k);
    fineSoln.initialize(f_foo, dx[1], k);
    fineData.setToZero();
    
    // COMPUTE INTERPOLATION
    crseData.copyTo(tempData);
    for (auto iter = fineLayout.begin(); iter.ok(); ++iter)
    {
        auto& tempData_i = tempData[*iter];
        auto& fineData_i = fineData[*iter];
        fineData_i |= interpolate(tempData_i);
    }

    // COMPUTE ERROR
    for (auto iter = fineLayout.begin(); iter.ok(); ++iter)
    {
        auto& fineData_i = fineData[*iter];
        auto& solnData_i = solnData[*iter];
        solnData_i -= fineData_i;
        solnData_i *= -1;
    }

    // WRITE OUTPUTS
    HDF5Handler h5;
    h5.writeLevel(dx[0], crseData, "CoarseData");    
    h5.writeLevel(dx[1], fineData, "FineData");    
    h5.writeLevel(dx[1], solnData, "Error");    
}
