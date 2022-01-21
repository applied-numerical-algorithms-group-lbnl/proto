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

int main(int argc, char** argv)
{
    // SETUP
    HDF5Handler h5;
    InputArgs args;
    args.parse();
    args.print();

    int domainSize = 64;
    int boxSize = 16;
    int numIter = 3;
    int ghostSize = 1;
    double k = 1;
    double s = 0.25;
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
        double dx = physDomainSize / domainSize;
        Point origin = Point::Ones(domainSize / 2) ;
        Box domainBox = Box::Cube(domainSize);

        BoxOp_Laplace<double> op(dx);
        
        BoxData<double> Phi_0(domainBox.grow(3));                
        BoxData<double> Phi(domainBox.grow(2));                
        BoxData<double> LPhi(domainBox);                
        BoxData<double> LPhiSln_0(domainBox.grow(1));                
        BoxData<double> LPhiSln(domainBox);                
        BoxData<double> LPhiErr(domainBox);                
        BoxData<double> RhoInv_0(domainBox.grow(3));                
        BoxData<double> RhoInv(domainBox.grow(2));                
        BoxData<double, DIM> RhoInv_f(domainBox.grow(1));                

        forallInPlace_p(f_Phi, Phi_0, dx, origin, k);
        Operator::convolve(Phi, Phi_0);
        forallInPlace_p(f_LRhoPhi, LPhiSln_0, dx, origin, k);
        Operator::convolve(LPhiSln, LPhiSln_0);
        forallInPlace_p(f_Rho, RhoInv_0, dx, origin, k);
        Operator::convolve(RhoInv, RhoInv_0);
        for (int dir = 0; dir < DIM; dir++)
        {
            auto interp = Stencil<double>::CellToFace(dir, Side::Lo);
            auto rho_d  = slice(RhoInv_f, dir);
            rho_d |= interp(RhoInv);
        }
        LPhi.setToZero();
        LPhiErr.setToZero();

        h5.writePatch(dx, Phi, "PHI_N%i", nn);
        h5.writePatch(dx, LPhiSln, "LPHI_SLN_N%i", nn);
        
        //op(LPhi, Phi, RhoInv_f);
        op(LPhi, Phi);
        h5.writePatch(dx, LPhi, "LPHI_N%i", nn);

        LPhi.copyTo(LPhiErr);
        LPhiErr -= LPhiSln;
        err[nn] = LPhiErr.absMax();
        h5.writePatch(dx, LPhiErr, "LPhi_ERR_N%i", nn);
        
        std::cout << "Error: " << err[nn] << std::endl;
        domainSize *= 2;
    }
        
    for (int ii = 1; ii < numIter; ii++)
    {
        std::cout << "Convergence Rate: " << log(err[ii-1] / err[ii]) / log(2.0) << std::endl;
    }

    return 0;
}

