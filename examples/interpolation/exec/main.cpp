#include "Proto.H"
#include "InputParser.H"

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
    InputArgs args;
    args.parse(); 
    args.print(); 
       
    // GENERAL PARAMETERS
    int numIter = 3;
    int domainSize = 64;
    int boxSize = 8;
    int refRatio = 2;
    int testNum = 0;
    int interpOrder = 2;
    
    args.set("numIter",         &numIter);
    args.set("domainSize",      &domainSize);
    args.set("boxSize",         &boxSize);
    args.set("refRatio",        &refRatio);
    args.set("testNum",         &testNum);
    args.set("interpOrder",     &interpOrder);
   
    
    // INTERP STENCIL PARAMETERS

    // Kernel: box which contains the entire footprint of the interp stencil 
    Box interpStencilShiftKernel = Box::Kernel(1);
    // MaxShift: defines how much of the kernel will be included in the stencil footprint.
    int interpStencilMaxShift = 2;
    // PolynomialOrder: Order of interpolant.
    int interpStencilPolynomialOrder = interpOrder;

    std::cout << "Test " << testNum << ": Interp Order = " << interpOrder << std::endl;

    switch (testNum)
    {
        case 0:
        {
#if DIM == 2
            std::cout << "Stencil Footprint (2D): " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "|   |   | 2 |   |   | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "|   | 2 | 1 | 2 |   | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "| 2 | 1 | 0 | 1 | 2 | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "|   | 2 | 1 | 2 |   | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "|   |   | 2 |   |   | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << std::endl;
            std::cout << "Numbers are shift-distance from origin" << std::endl;
            std::cout << "Empty cells are not part of the stencil footprint" << std::endl;
#endif
            interpStencilMaxShift = 2;
            interpStencilShiftKernel = Box::Kernel(2);
            break;
        }
        case 1:
        {
#if DIM == 2
            std::cout << "Stencil Footprint (2D): " << std::endl;
            std::cout << "+---+---+---+---+ " << std::endl;
            std::cout << "|   | 2 |   |   | " << std::endl;
            std::cout << "+---+---+---+---+ " << std::endl;
            std::cout << "| 2 | 1 | 1 |   | " << std::endl;
            std::cout << "+---+---+---+---+ " << std::endl;
            std::cout << "| 1 | 0 | 1 | 2 | " << std::endl;
            std::cout << "+---+---+---+---+ " << std::endl;
            std::cout << "| 2 | 1 | 1 |   | " << std::endl;
            std::cout << "+---+---+---+---+ " << std::endl;
            std::cout << "|   | 2 |   |   | " << std::endl;
            std::cout << "+---+---+---+---+ " << std::endl;
            std::cout << "Numbers are shift-distance from origin" << std::endl;
            std::cout << "Empty cells are not part of the stencil footprint" << std::endl;
#endif
        
            interpStencilMaxShift = 2;
            interpStencilShiftKernel = Box(Point(-1,-2), Point(2, 2));
            break;
        }
        case 2:
        {
#if DIM == 2
            std::cout << "Stencil Footprint (2D): " << std::endl;
            std::cout << "+---+---+---+ " << std::endl;
            std::cout << "| 2 |   |   | " << std::endl;
            std::cout << "+---+---+---+ " << std::endl;
            std::cout << "| 1 | 2 |   | " << std::endl;
            std::cout << "+---+---+---+ " << std::endl;
            std::cout << "| 0 | 1 | 2 | " << std::endl;
            std::cout << "+---+---+---+ " << std::endl;
            std::cout << "| 1 | 2 |   | " << std::endl;
            std::cout << "+---+---+---+ " << std::endl;
            std::cout << "| 2 |   |   | " << std::endl;
            std::cout << "+---+---+---+ " << std::endl;
            std::cout << "Numbers are shift-distance from origin" << std::endl;
            std::cout << "Empty cells are not part of the stencil footprint" << std::endl;
#endif
        
            interpStencilMaxShift = 2;
            interpStencilShiftKernel = Box(Point(0,-2), Point(2, 2));
            break;
        }
        case 3:
        {
#if DIM == 2
            std::cout << "Stencil Footprint (2D): " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "|   | 3 | 2 | 3 |   | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "| 3 | 2 | 1 | 2 | 3 | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "| 2 | 1 | 0 | 1 | 2 | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "| 3 | 2 | 1 | 2 | 3 | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "|   | 3 | 2 | 3 |   | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "Numbers are shift-distance from origin" << std::endl;
            std::cout << "Empty cells are not part of the stencil footprint" << std::endl;
#endif
        
            interpStencilMaxShift = 3;
            interpStencilShiftKernel = Box::Kernel(2);
            break;
        }
        case 4:
        {
#if DIM == 2
            std::cout << "Stencil Footprint (2D): " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "| 4 | 3 | 2 | 3 | 4 | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "| 3 | 2 | 1 | 2 | 3 | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "| 2 | 1 | 0 | 1 | 2 | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "| 3 | 2 | 1 | 2 | 3 | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "| 4 | 3 | 2 | 3 | 4 | " << std::endl;
            std::cout << "+---+---+---+---+---+ " << std::endl;
            std::cout << "Numbers are shift-distance from origin" << std::endl;
            std::cout << "Empty cells are not part of the stencil footprint" << std::endl;
#endif
        
            interpStencilMaxShift = 4;
            interpStencilShiftKernel = Box::Kernel(2);
            break;
        }
    }
    // BUILD INTERPOLATION OPERATOR
    std::cout << std::endl;
    std::cout << "Interpolating Polynomial Order: " << interpOrder << std::endl;
    std::cout << "InterpStencil Footprint Kernel: " << interpStencilShiftKernel << std::endl;
    std::cout << "InterpStencil Max Shift: " << interpStencilMaxShift << std::endl;    
    auto interpolate = InterpStencil<double>::Build(
            interpStencilMaxShift,
            interpStencilShiftKernel,
            interpStencilPolynomialOrder,
            refRatio);
   
    std::cout << "Span of interpolation stencil: " << interpolate.ghost() << std::endl; 
    double physDomainSize = 1.0;
    double k = 1.0;

    double error[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx[DIM];
        dx[0] = physDomainSize / domainSize;
        dx[1] = dx[0] / refRatio;


        // BUILD GEOMETRY
        Point boxSizeVect = Point::Ones(boxSize);
        Point refRatioVect = Point::Ones(refRatio);
        Box domainBox = Box::Cube(domainSize);
        std::array<bool, DIM> periodicity;
        for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }
        Box refinedRegion = Box::Cube(domainSize / 2).shift(Point::Ones(domainSize / 4)).refine(refRatio);
        Box refinedPatches = refinedRegion.coarsen(boxSize);
        std::vector<Point> refinedPatchPoints;
        for (auto iter = refinedPatches.begin(); iter.ok(); ++iter)
        {
            refinedPatchPoints.push_back(*iter);
        }

        ProblemDomain crseDomain(domainBox, periodicity);
        ProblemDomain fineDomain = crseDomain.refine(refRatioVect);
        DisjointBoxLayout crseLayout(crseDomain, boxSizeVect);
        DisjointBoxLayout fineLayout(fineDomain, refinedPatchPoints, boxSizeVect);
        DisjointBoxLayout crseFineLayout = fineLayout.coarsen(refRatioVect);    

        // INITIALIZE DATA HOLDERS
        LevelBoxData<double> crseData(crseLayout, Point::Ones());
        LevelBoxData<double> tempData(crseFineLayout, Point::Ones() + interpolate.ghost());
        LevelBoxData<double> fineData(fineLayout, Point::Ones());
        LevelBoxData<double> fineSoln(fineLayout, Point::Ones());

        crseData.initConvolve(f_foo, dx[0], k);
        fineSoln.initConvolve(f_foo, dx[1], k);
        fineData.setToZero();

        // COMPUTE INTERPOLATION

        crseData.copyTo(tempData);
        //      The temporary is needed to guarantee 1:1 correspondence between coarse and fine patches
        //      when applying the stencil in the following loop
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
            auto& fineSoln_i = fineSoln[*iter];
            fineSoln_i -= fineData_i;
            fineSoln_i *= -1;
        }
        error[nn] = fineSoln.absMax();
        std::cout << "Error: " << error[nn] << std::endl;

        // WRITE OUTPUTS
        HDF5Handler h5;
        h5.writeLevel(dx[0], tempData, "TempData_D%i", domainSize);
        h5.writeLevel(dx[0], crseData, "CoarseData_D%i", domainSize);    
        h5.writeLevel(dx[1], fineData, "FineData_D%i", domainSize);    
        h5.writeLevel(dx[1], fineSoln, "Error_D%i", domainSize);    
        
        domainSize *= 2;
    }

    for (int ii = 1; ii < numIter; ii++)
    {
        std::cout << "Convergence Rate: " << log(error[ii-1]/error[ii])/log(2.0) << std::endl;
    }
}
