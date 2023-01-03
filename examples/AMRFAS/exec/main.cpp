#include "Proto.H"
#include "InputParser.H"
#include "BoxOp_Laplace.H"
#include "AMRSolver_FASMultigrid.H"

using namespace Proto;

PROTO_KERNEL_START void f_force_avg_0(const Point& a_pt, Var<double> a_data, Array<double, DIM> a_dx)
{
    double x0[DIM];
    double x1[DIM];
    double* dx = a_dx.data();
    double a = 0.125;
    for (int dir = 0; dir < DIM; dir++)
    {
        x0[dir] = a_pt[dir]*dx[dir] + a;
        x1[dir] = x0[dir] + dx[dir];
    }
    
    double k = M_PI*2;
    a_data(0) = + cos(k*x1[0])*cos(k*x1[1])
                - cos(k*x0[0])*cos(k*x1[1])
                - cos(k*x1[0])*cos(k*x0[1])
                + cos(k*x0[0])*cos(k*x0[1]);
    a_data(0) *= 1.0/(k*k*dx[0]*dx[1]);
}
PROTO_KERNEL_END(f_force_avg_0, f_force_avg);

PROTO_KERNEL_START void f_soln_avg_0(const Point& a_pt, Var<double> a_data, Array<double, DIM> a_dx)
{
    double x0[DIM];
    double x1[DIM];
    double* dx = a_dx.data();
    double a = 0.125;
    for (int dir = 0; dir < DIM; dir++)
    {
        x0[dir] = a_pt[dir]*dx[dir] + a;
        x1[dir] = x0[dir] + dx[dir];
    }
    
    double k = M_PI*2;
    a_data(0) = + cos(k*x1[0])*cos(k*x1[1])
                - cos(k*x0[0])*cos(k*x1[1])
                - cos(k*x1[0])*cos(k*x0[1])
                + cos(k*x0[0])*cos(k*x0[1]);
    a_data(0) *= 1.0/(k*k*dx[0]*dx[1]);
    a_data(0) *= -1.0/(DIM*pow(k, 2.0));
}
PROTO_KERNEL_END(f_soln_avg_0, f_soln_avg);

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif
   
    // SETUP
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    using Proto::pout;
    typedef BoxOp_Laplace<double> OP;

    int domainSize = 64;
    int boxSize = 32;
    int numIter = 2;
    int numLevels = 2;
    int solveIter = 20;
    double tolerance = 1e-10;
    int refRatio = 4;
    Array<bool, DIM> periodicity;
    periodicity.fill(true);
    
    InputArgs args;
    args.add("domainSize", domainSize);
    args.add("boxSize",    boxSize);
    args.add("numIter",    numIter);
    args.add("numLevels",  numLevels);
    args.add("solveIter",  solveIter);
    args.add("tolerance",  tolerance);
    args.add("refRatio",   refRatio);
    args.add("periodic_x", periodicity[0]);
#if DIM > 1
    args.add("periodic_y", periodicity[1]);
#endif
#if DIM > 2
    args.add("periodic_z", periodicity[2]);
#endif
    args.parse(argc, argv);
    args.print();
    PR_TIMER_SETFILE(to_string(domainSize) + ".DIM" + to_string(DIM) + ".AMRFAS.time.table");
    double physDomainSize = 1.0;

    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        // GRIDS
        Array<double, DIM> dx;
        dx.fill(physDomainSize / domainSize);
        
        std::vector<Point> refRatios;
        refRatios.resize(numLevels-1, Point::Ones(refRatio));
        std::vector<DisjointBoxLayout> layouts;
        layouts.resize(numLevels);

        Point boxSizeVect = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain domain(domainBox, periodicity);
        layouts[0].define(domain, boxSizeVect);

        Box refinedRegion = domainBox;
        for (int lvl = 1; lvl < numLevels; lvl++)
        {
            refinedRegion = refinedRegion.grow(-refinedRegion.sizes() / 4).refine(refRatios[lvl-1]);
            ProblemDomain fineDomain = layouts[lvl-1].domain().refine(refRatios[lvl-1]);
            layouts[lvl].define(fineDomain, refinedRegion, boxSizeVect);
        }
        AMRGrid grid(layouts, refRatios, numLevels);

        // SOLVER
        AMRSolver_FASMultigrid<BoxOp_Laplace, double> solver(grid, dx);

        // DATA
        AMRData<double> Phi(grid,    OP::ghost()); 
        AMRData<double> G(grid,      Point::Zeros());
        AMRData<double> PhiSln(grid, OP::ghost()); 
        AMRData<double> PhiErr(grid, Point::Zeros()); 
        AMRData<double> Res(grid,    Point::Zeros());

        Phi.setToZero();
        Array<double,DIM> dev_x(dx);
        G.initialize(dev_x, f_force_avg);
        PhiSln.initialize(dev_x, f_soln_avg);

#ifdef PR_HDF5
        h5.writeAMRData(dx, PhiSln, "SLN_N%i", nn);
        h5.writeAMRData(dx, G, "RHS_N%i", nn);
#endif
        
        // SOLVE
        pout() << "Integral of RHS: " << G.integrate(dx) << std::endl;
        {
          PR_TIMERS("AMR FAS solve");
          solver.solve(Phi, G, solveIter, tolerance);
        }
        Phi.averageDown();
        double phiav = Phi.integrate(dx)/pow(physDomainSize,DIM);
        Phi += (-phiav);
        if (procID() == 0)
        {
            std::cout << "average of phi: " << phiav << std::endl;
        }
        // COMPUTE ERROR
        PhiErr.setToZero();
        Phi.copyTo(PhiErr);
        PhiErr.increment(PhiSln, -1);
        PhiErr.averageDown();
        
        //err[nn] = PhiErr.integrateAbs(dx);
        err[nn] = PhiErr.absMax();
        //pout() << "Error: " << err[nn] << std::endl;
        if (procID() == 0)
        {
            std::cout << "Error: " << err[nn] << std::endl;
        }

        domainSize *= 2;
    }
        
    for (int ii = 1; ii < numIter; ii++)
    {
      //pout() << "Convergence Rate: " << log(err[ii-1] / err[ii]) / log(2.0) << std::endl;
        if (procID() == 0)
        {
            std::cout << "Convergence Rate: " << log(err[ii-1] / err[ii]) / log(2.0) << std::endl;
        }
    }
    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
    #endif
}
