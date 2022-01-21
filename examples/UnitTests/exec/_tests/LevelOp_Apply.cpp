#include "Proto.H"
#include "InputParser.H"
#include "BoxOp_Laplace.H"
#include "func.H"

using namespace Proto;

void f_LRhoPhi(Point& a_pt, Var<double>& a_data, double a_dx,
    Point a_origin, double a_k)
{
    double x[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = (a_pt[ii] - a_origin[ii])*a_dx + a_dx / 2.0;
    }
    double k1 = 2.0*M_PI*a_k;
    double k2 = 4.0*M_PI*a_k;
    double s1 = sin(k1*(x[0] + x[1]));
    double s2 = sin(k2*x[0]);
    double c1 = cos(k1*(x[0] + x[1]));
    double c2 = cos(k2*x[0]);
   
    //a_data(0) = k2*c2*k1*c1 - 2*s2*k1*k1*s1;
    a_data(0) = -DIM*pow(k1, 2)*s1;
}

void f_Phi(Point& a_pt, Var<double>& a_data, double a_dx,
    Point a_origin, double a_k)
{
    double x[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = (a_pt[ii] - a_origin[ii])*a_dx + a_dx / 2.0;
    }
    
    a_data(0) = sin(2.0*M_PI*a_k*(x[0] + x[1]));
}

void f_Rho(Point& a_pt, Var<double>& a_data, double a_dx,
    Point a_origin, double a_k)
{
    double x[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = (a_pt[ii] - a_origin[ii])*a_dx + a_dx / 2.0;
    }
    
    a_data(0) = sin(4.0*M_PI*a_k*(x[0]));

}

void f_LPhi(Point& a_pt, Var<double>& a_data, double a_dx,
    Point a_origin, double a_k)
{

    double x[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = (a_pt[ii] - a_origin[ii])*a_dx + a_dx / 2.0;
    }
    double k = 2*M_PI*a_k;
    double s = sin(k*(x[0] + x[1]));
    
    a_data(0) = -2*pow(k, 2)*s;
}

int main(int argc, char** argv)
{
    // SETUP
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif
    
    HDF5Handler h5;
    InputArgs args;
    args.parse();
    args.print();

    int domainSize = 64;
    int boxSize = 16;
    int numIter = 3;
    int ghostSize = 1;
    double k = 1;
    double s = 0.1;
    double physDomainSize = 1;
    int refRatio = PR_AMR_REFRATIO;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    args.set("domainSize", &domainSize);
    args.set("boxSize",    &boxSize);
    args.set("numIter",    &numIter);
    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);

    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        // BUILD GRIDS
        double dx = physDomainSize / domainSize;
        Point origin = Point::Ones(domainSize / 2.0);
        
        Point boxSizeV = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain domain(domainBox, periodicity);
        DisjointBoxLayout layout(domain, boxSizeV);
        
        LevelBoxData<double>      Phi(layout,      Point::Ones(3));
        LevelBoxData<double>      LPhi(layout,     Point::Zeros());
        LevelBoxData<double>      LPhiSln(layout,  Point::Zeros());
        LevelBoxData<double>      LPhiErr(layout,  Point::Zeros());
        LevelBoxData<double>      RhoInv(layout,   Point::Ones(3));
        LevelBoxData<double, DIM> RhoInv_f(layout, Point::Ones(2));
        
        Phi.initConvolve(f_Phi, dx, origin, k);
        LPhiSln.initConvolve(f_LRhoPhi, dx, origin, k);
        RhoInv.initConvolve(f_Rho, dx, origin, k);
        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            auto& rho_i = RhoInv[*iter];
            auto& rho_f_i = RhoInv_f[*iter];
            for (int dir = 0; dir < DIM; dir++)
            {
                auto interp = Stencil<double>::CellToFace(dir, Side::Lo);
                auto rho_f_id = slice(rho_f_i, dir);
                rho_f_id |= interp(rho_i);
            }
        }
        RhoInv_f.exchange();
        LPhi.setToZero();
        LPhiErr.setToZero();
        
        h5.writeLevel(dx, Phi,      "PHI");
        h5.writeLevel(dx, LPhiSln,  "LPHI_SLN");
        h5.writeLevel(dx, RhoInv,   "RHO");
        h5.writeLevel(dx, RhoInv_f, "RHO_F");
       
        LevelOp<BoxOp_Laplace, double> op(dx);
        //op(LPhi, Phi, RhoInv_f);
        op(LPhi, Phi);

        h5.writeLevel(dx, LPhi, "LPHI");
        
        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            auto& lphi_i = LPhi[*iter];
            auto& sln_i = LPhiSln[*iter];
            auto& err_i = LPhiErr[*iter];
            
            lphi_i.copyTo(err_i);
            err_i -= sln_i;
        }
        err[nn] = LPhiErr.absMax();
        h5.writeLevel(dx, LPhiErr, "LPHI_ERR");

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

