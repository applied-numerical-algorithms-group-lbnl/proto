#include "Proto.H"
#include "InputParser.H"
#include "BoxOp_Laplace.H"
#include "LevelSolver_FASMultigrid.H"
#include "func.H"

using namespace Proto;

PROTO_KERNEL_START void f_force_0(const Point& a_pt, Var<double> a_data, double a_dx)
{
    double x[DIM];
    
    a_data(0) = 1.0;
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = a_pt[dir]*a_dx + a_dx/2.0;
        a_data(0) = a_data(0)*sin(M_PI*2*(x[dir] + .125));
    }
}
PROTO_KERNEL_END(f_force_0, f_force);

PROTO_KERNEL_START void f_soln_0(const Point& a_pt, Var<double> a_data, double a_dx)
{
    double x[DIM];
    
    a_data(0) = -1/(DIM*pow(2.0*M_PI, 2));
    for (int dir = 0; dir < DIM; dir++)
    {
        x[dir] = a_pt[dir]*a_dx + a_dx/2.0;
        a_data(0) = a_data(0)*sin(M_PI*2*(x[dir] + .125));
    }
}
PROTO_KERNEL_END(f_soln_0, f_soln);

int MG_LEVEL = 0;
int SOLVE_ITER = 0;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    HDF5Handler h5;
    InputArgs args;
    args.parse();
    args.print();

    int domainSize = 64;
    int boxSize = 16;
    int numIter = 3;
    int solveIter = 20;
    int numLevels = log(1.0*domainSize)/log(2.0) + 1;
    double tolerance = 1e-10;
    int ghostSize = 1;
    double k = 1;
    double physDomainSize = 1;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    args.set("domainSize", &domainSize);
    args.set("boxSize",    &boxSize);
    args.set("numIter",    &numIter);
    args.set("numLevels",  &numLevels);
    args.set("solveIter",  &solveIter);
    args.set("tolerance",  &tolerance);
    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);
    
    typedef BoxOp_Laplace<double> OP;
    
    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        // BUILD GRIDS
        double dx = physDomainSize / domainSize;
        // coarse layout
        Point boxSizeV = Point::Ones(boxSize);
        Box crseDomainBox = Box::Cube(domainSize);
        ProblemDomain crseDomain(crseDomainBox, periodicity);
        DisjointBoxLayout crseLayout(crseDomain, boxSizeV);
        // fine layout 
        int amrRefRatio = pow(PR_MG_REFRATIO, numLevels);
        Box fineDomainBox = Box::Cube(domainSize / 2).shift(Point::Ones(domainSize / 4));
        fineDomainBox = fineDomainBox.refine(amrRefRatio);
        ProblemDomain fineDomain = crseDomain.refine(Point::Ones(amrRefRatio));
        DisjointBoxLayout fineLayout(fineDomain, fineDomainBox, boxSizeV);
    
        
        double cdx = dx;
        double fdx = dx / amrRefRatio;
        // solver
        LevelSolver_FASMultigrid<BoxOp_Laplace, double> solver(fineLayout, numLevels, fdx);

        // data holders
        LevelBoxData<double, OP::numState()> Phi(fineLayout, OP::ghost());
        LevelBoxData<double, OP::numState()> PhiSln(fineLayout, Point::Zeros());
        LevelBoxData<double, OP::numState()> PhiErr(fineLayout, Point::Zeros());
        LevelBoxData<double, OP::numState()> G(fineLayout, Point::Zeros());
        LevelBoxData<double, OP::numState()> PhiBC(crseLayout, OP::ghost());

        Phi.setToZero();
        G.initConvolve(f_force, fdx);
        PhiSln.initConvolve(f_soln, fdx);
        PhiBC.initConvolve(f_soln, cdx);
        PhiErr.setToZero();
       
        h5.writeLevel(fdx, PhiSln, "SLN");
        h5.writeLevel(cdx, PhiBC,  "PHI_BC");
        h5.writeLevel(fdx, Phi,     "PHI_0");
        h5.writeLevel(fdx, G,       "G");

        solver.defineAsSubcycler(fineLayout, numLevels, fdx);
        solver.interpBC(PhiBC);
        solver.solve(Phi, G, solveIter, tolerance);

        h5.writeLevel(dx, Phi,    "PHI_1");
        
        for (auto iter = fineLayout.begin(); iter.ok(); ++iter)
        {
            auto& phi_i = Phi[*iter];
            auto& sln_i = PhiSln[*iter];
            auto& err_i = PhiErr[*iter];
        
            phi_i.copyTo(err_i);
            err_i -= sln_i;
        }
        h5.writeLevel(dx, PhiErr,    "ERR");
        err[nn] = PhiErr.absMax();
        pout() << "Error: " << err[nn] << std::endl;
        domainSize *= 2;
    }
        
    for (int ii = 1; ii < numIter; ii++)
    {
        pout() << "Convergence Rate: " << log(err[ii-1] / err[ii]) / log(2.0) << std::endl;
    }

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

