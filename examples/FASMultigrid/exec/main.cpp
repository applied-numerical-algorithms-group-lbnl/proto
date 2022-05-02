#include "ProtoAMR.H"
#include "InputParser.H"
#include "BoxOp_Laplace.H"
#include "LevelSolver_FASMultigrid.H"

using namespace Proto;

PROTO_KERNEL_START
void f_force_avg_0(const Point& a_pt, Var<double> a_data, double a_dx)
{
    double x0[DIM];
    double x1[DIM];
    
    double a = 0.125;
    for (int dir = 0; dir < DIM; dir++)
    {
        x0[dir] = a_pt[dir]*a_dx + a;
        x1[dir] = x0[dir] + a_dx;
    }
    
    double k = M_PI*2;
    a_data(0) = + cos(k*x1[0])*cos(k*x1[1])
                - cos(k*x0[0])*cos(k*x1[1])
                - cos(k*x1[0])*cos(k*x0[1])
                + cos(k*x0[0])*cos(k*x0[1]);
    a_data(0) *= 1.0/(k*k*a_dx*a_dx);
}
PROTO_KERNEL_END(f_force_avg_0, f_force_avg);

PROTO_KERNEL_START
void f_soln_avg_0(const Point& a_pt, Var<double> a_data, double a_dx)
{
    double x0[DIM];
    double x1[DIM];
    
    double a = 0.125;
    for (int dir = 0; dir < DIM; dir++)
    {
        x0[dir] = a_pt[dir]*a_dx + a;
        x1[dir] = x0[dir] + a_dx;
    }
    
    double k = M_PI*2;
    a_data(0) = + cos(k*x1[0])*cos(k*x1[1])
                - cos(k*x0[0])*cos(k*x1[1])
                - cos(k*x1[0])*cos(k*x0[1])
                + cos(k*x0[0])*cos(k*x0[1]);
    a_data(0) *= 1.0/(k*k*a_dx*a_dx);
    a_data(0) *= -1.0/(DIM*pow(k, 2.0));
}
PROTO_KERNEL_END(f_soln_avg_0, f_soln_avg);

int MG_LEVEL = 0;
int SOLVE_ITER = 0;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    using Proto::pout;
    HDF5Handler h5;

    int domainSize = 64;
    int boxSize = 16;
    int numIter = 3;
    int numLevels = log(1.0*domainSize)/log(2.0) + 1;
    int solveIter = 20;
    double tolerance = 1e-10;
    int refRatio = 2;
    std::array<bool, DIM> periodicity;
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
    args.add("periodic_y", periodicity[1]);
    args.parse(argc, argv);
    args.print();
    
    double physDomainSize = 1;
    
    typedef BoxOp_Laplace<double> OP;
    
    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        // BUILD GRIDS
        
        double dx = physDomainSize / domainSize;
        std::array<double, DIM> dxVect;
        dxVect.fill(dx);

        // coarse layout
        Point boxSizeV = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain domain(domainBox, periodicity);
        DisjointBoxLayout layout(domain, boxSizeV);
        
        // solver
        LevelSolver_FASMultigrid<BoxOp_Laplace, double> solver(layout, Point::Ones(refRatio), numLevels, dx);
        solver.setVerbose(true);

        // data holders
        LevelBoxData<double, OP::numState()> Phi(layout, OP::ghost());
        LevelBoxData<double, OP::numState()> PhiSln(layout, Point::Zeros());
        LevelBoxData<double, OP::numState()> PhiErr(layout, Point::Zeros());
        LevelBoxData<double, OP::numState()> G(layout, Point::Zeros());
        Phi.setToZero();
        G.initialize(f_force_avg, dx);
        PhiSln.initialize(f_soln_avg, dx);
       
        h5.writeLevel(dx, PhiSln, "SLN_N%i", nn);
        h5.writeLevel(dx, G,       "G_N%i", nn);
        pout() << "Integral of RHS: " << G.integrate(dx) << std::endl;

        solver.solve(Phi, G, solveIter, tolerance);

        h5.writeLevel(dx, Phi,    "PHI_N%i", nn);
        
        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            auto& phi_i = Phi[*iter];
            auto& sln_i = PhiSln[*iter];
            auto& err_i = PhiErr[*iter];
        
            phi_i.copyTo(err_i);
            err_i -= sln_i;
        }
        h5.writeLevel(dx, PhiErr,    "ERR_N%i", nn);
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

