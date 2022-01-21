#include "Proto.H"
#include "InputParser.H"
#include "BoxOp_Laplace.H"
#include "func.H"

using namespace Proto;

void f_LRhoPhi(Point& a_pt, Var<double>& a_data, double a_dx, double a_k)
{
    double x[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = a_pt[ii]*a_dx + a_dx / 2.0;
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

void f_Phi(Point& a_pt, Var<double>& a_data, double a_dx, double a_k)
{
    double x[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = a_pt[ii]*a_dx + a_dx / 2.0;
    }
    
    a_data(0) = sin(2.0*M_PI*a_k*(x[0] + x[1]));
}

void f_Rho(Point& a_pt, Var<double>& a_data, double a_dx, double a_k)
{
    double x[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = a_pt[ii]*a_dx + a_dx / 2.0;
    }
    
    a_data(0) = sin(4.0*M_PI*a_k*(x[0]));

}

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
    int numLevels = 2;
    int ghostSize = 1;
    double k = 1;
    double physDomainSize = 1;
    int refRatio = PR_AMR_REFRATIO;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    args.set("domainSize", &domainSize);
    args.set("boxSize",    &boxSize);
    args.set("numIter",    &numIter);
    args.set("numLevels",  &numLevels);
    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);
    
    double err[numIter];
    double err0[numIter];
    double err1[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        // BUILD GRIDS
        double dx = physDomainSize / domainSize;
        
        std::vector<DisjointBoxLayout> layouts;
        layouts.resize(numLevels);
                
        // coarse layout
        Point boxSizeV = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain domain(domainBox, periodicity);
        layouts[0].define(domain, boxSizeV);

        Box refinedRegion = domainBox;
        for (int lvl = 1; lvl < numLevels; lvl++)
        {
            Point prevSize = refinedRegion.high() - refinedRegion.low() + Point::Ones();
            refinedRegion = refinedRegion.grow(-prevSize / 4).refine(refRatio);
            
            Box refinedRegionPatches = refinedRegion.coarsen(boxSizeV);
            std::vector<Point> fineLayoutPatches;
            bool skip = true;
            for (auto iter = refinedRegionPatches.begin(); iter.ok(); ++iter)
            {
                if (skip)
                {
                    skip = false;
                    continue;
                }
                fineLayoutPatches.push_back(*iter);
            }
            ProblemDomain fineDomain = layouts[lvl-1].domain().refine(Point::Ones(refRatio));
            layouts[lvl].define(fineDomain, fineLayoutPatches, boxSizeV);
        }

        AMRGrid grid(layouts, numLevels);
        
        // TEST CODE
        AMRData<double>      Phi(grid,      Point::Ones(3));
        AMRData<double>      LPhi(grid,     Point::Zeros());
        AMRData<double>      LPhiSln(grid,  Point::Zeros());
        AMRData<double>      LPhiErr(grid,  Point::Zeros());
        AMRData<double>      RhoInv(grid,   Point::Ones(3));
        AMRData<double, DIM> RhoInv_f(grid, Point::Ones(2));
        LevelBoxData<double> LPhi0(grid[0], Point::Zeros());
        LevelBoxData<double> LPhiErr0(grid[0], Point::Zeros());
        LevelBoxData<double> LPhi1(grid[1], Point::Zeros());
        LevelBoxData<double> LPhiErr1(grid[1], Point::Zeros());

        Phi.initConvolve(dx, f_Phi, k);
        LPhiSln.initConvolve(dx, f_LRhoPhi, k);
        RhoInv.initConvolve(dx, f_Rho, k);
        
        double dx_lvl = dx;
        for (int lvl = 0; lvl < numLevels; lvl++)
        {
            for (auto iter = grid[lvl].begin(); iter.ok(); ++iter)
            {
                auto& rho_i = RhoInv[lvl][*iter];
                auto& rho_f_i = RhoInv_f[lvl][*iter];
                for (int dir = 0; dir < DIM; dir++)
                {
                    auto interp = Stencil<double>::CellToFace(dir, Side::Lo);
                    auto rho_f_id = slice(rho_f_i, dir);
                    rho_f_id |= interp(rho_i);
                }
            }
        }
        RhoInv_f.exchange();
        LPhi.setToZero();
        LPhi0.setToZero();
        LPhi1.setToZero();
        LPhiErr.setToZero();
        LPhiErr0.setToZero(); 
        LPhiErr1.setToZero(); 

        h5.writeAMRData(dx, Phi,      "PHI_N%i", nn);
        h5.writeAMRData(dx, LPhiSln,  "LPHI_SLN_N%i", nn);
        h5.writeAMRData(dx, RhoInv,   "RHO");
        h5.writeAMRData(dx, RhoInv_f, "RHO_F");
      
        Phi.averageDown();
        auto INTERP = InterpStencil<double>::Build(4, Box::Kernel(2), 4, refRatio);
        interpBoundaries(Phi[0], Phi[1], INTERP);
        Phi.exchange();
        AMROp<BoxOp_Laplace, double> op(grid, dx);
        //op(LPhi, Phi, RhoInv_f);
        op(LPhi, Phi);
        op.levelApply(LPhi0, Phi, 0);
        op.levelApply(LPhi1, Phi, 1);
        LPhi.averageDown(); 
          
        h5.writeAMRData(dx, LPhi, "LPHI");
        h5.writeLevel(dx, LPhi0, "LPHI0");
        
        for (int lvl = 0; lvl < numLevels; lvl++)
        {
            for (auto iter = grid[lvl].begin(); iter.ok(); ++iter)
            {
                auto& lphi_i = LPhi[lvl][*iter];
                auto& sln_i = LPhiSln[lvl][*iter];
                auto& err_i = LPhiErr[lvl][*iter];

                lphi_i.copyTo(err_i);
                err_i -= sln_i;
            }
        }

        for (auto iter = grid[0].begin(); iter.ok(); ++iter)
        {
            auto& lphi_i = LPhi0[*iter];
            auto& sln_i = LPhiSln[0][*iter];
            auto& err_i = LPhiErr0[*iter];

            lphi_i.copyTo(err_i);
            err_i -= sln_i;
        }
        for (auto iter = grid[1].begin(); iter.ok(); ++iter)
        {
            auto& lphi_i = LPhi1[*iter];
            auto& sln_i = LPhiSln[1][*iter];
            auto& err_i = LPhiErr1[*iter];

            lphi_i.copyTo(err_i);
            err_i -= sln_i;
        }
        LPhiErr.averageDown();
        err[nn] = LPhiErr.integrateAbs(dx);
        err0[nn] = LPhiErr0.integrateAbs(dx);
        err1[nn] = LPhiErr1.integrateAbs(dx);
        h5.writeAMRData(dx, LPhiErr, "LPHI_ERR_N%i", nn);
        h5.writeLevel(dx, LPhiErr0, "LPHI0_ERR_N%i", nn);
        h5.writeLevel(dx / refRatio, LPhiErr1, "LPHI1_ERR_N%i", nn);
        
        pout() << "Error (AMR Apply):   " << err[nn] << std::endl;
        pout() << "Error (Level): L0: " << err0[nn] << " | L1: " << err1[nn] << std::endl;
        
        domainSize *= 2;
    }
        
    for (int ii = 1; ii < numIter; ii++)
    {
        pout() << "Convergence Rate (AMR Apply):   " << log(err[ii-1] / err[ii]) / log(2.0) << std::endl;
        pout() << "Convergence Rate (Level 0): " << log(err0[ii-1] / err0[ii]) / log(2.0) << std::endl;
        pout() << "Convergence Rate (Level 1): " << log(err1[ii-1] / err1[ii]) / log(2.0) << std::endl;
    }

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

