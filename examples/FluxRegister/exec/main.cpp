#include "Proto.H"
#include "InputParser.H"
#include <iomanip>

using namespace Proto;

void f_phi(Point& a_pt, Var<double>& a_phi, double a_dx, double a_kx, double a_ky)
{
    double x[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = a_pt[ii]*a_dx + a_dx / 2.0;
    }

    a_phi(0) = sin(4.0*M_PI*(a_kx*x[0] + a_ky*x[1]));
}

void f_lphi(Point& a_pt, Var<double>& a_phi, double a_dx, double a_kx, double a_ky)
{
    double x[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = a_pt[ii]*a_dx + a_dx / 2.0;
    }
      
    a_phi(0) = - pow(4.0*M_PI*a_kx, 2)*sin(4.0*M_PI*(a_kx*x[0] + a_ky*x[1]))
               - pow(4.0*M_PI*a_ky, 2)*sin(4.0*M_PI*(a_kx*x[0] + a_ky*x[1]));
}


int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init (&argc, &argv);
    setPoutBaseName("FluxRegister");
#endif

    // READ INPUT PARAMETERS
    HDF5Handler h5;
    
    int domainSize = 64;
    int boxSize = 32;
    int numIter = 3;
    double kx = 1;
    double ky = 1;

    InputArgs args;
    args.parse(); //assumes a file "inputs" exists.
    args.set("domainSize", &domainSize);
    args.set("boxSize", &boxSize);
    args.set("numIter", &numIter);
    args.set("kx", &kx);
    args.set("ky", &ky);
    
    int refRatio = PR_AMR_REFRATIO;
    PR_TIMER_SETFILE(to_string(domainSize) + ".proto.time.table");
    PR_TIMERS("main");

    double kk = sqrt(kx*kx + ky*ky);
    kx /= kk;
    ky /= kk;

    Point refRatioVect = Point::Ones(refRatio);
    Point boxSizeVect = Point::Ones(boxSize);
    
    double err[numIter];
    double err_cf[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        double physDomainSize = 1.0;
        double dx[2];
        dx[0] = physDomainSize / domainSize;
        dx[1] = dx[0] / refRatio;

        // BUILD LAYOUTS
        std::array<bool, DIM> periodicity;
        for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

        Box crseDomainBox = Box::Cube(domainSize);
        Box cfDomainBox = Box::Cube(domainSize/2).shift(Point::Ones(domainSize/4));
        Box fineDomainBox = cfDomainBox.refine(refRatio);

        ProblemDomain crseDomain(crseDomainBox, periodicity);
        ProblemDomain fineDomain = crseDomain.refine(refRatioVect);

        Box fineIgnoreBox = Box(Point::Zeros(), Point(domainSize / 2 - 1, domainSize - 1));
        fineIgnoreBox = fineIgnoreBox.refine(refRatio).coarsen(boxSizeVect);
        std::vector<Point> finePatchPoints;
        Box finePatchBox = fineDomainBox.coarsen(boxSizeVect);
        for (auto iter = finePatchBox.begin(); iter.ok(); ++iter)
        {
            if (fineIgnoreBox.contains(*iter)) { continue; }
            finePatchPoints.push_back(*iter);
        }
        std::vector<DisjointBoxLayout> grids;
        grids.resize(2);
        grids[0].define(crseDomain, boxSizeVect);
        grids[1].define(fineDomain, finePatchPoints, boxSizeVect);
        AMRGrid grid(grids, 2);
        auto& crseLayout = grid[0];
        auto& fineLayout = grid[1];
        DisjointBoxLayout cfLayout = fineLayout.coarsen(Point::Ones(refRatio));

        // BUILD DATA HOLDERS
        AMRData<double>         Phi(grid, Point::Ones(3));
        AMRData<double>         LPhi(grid, Point::Zeros());
        AMRData<double>         LPhiSoln(grid, Point::Zeros());
        AMRData<double, DIM>    Flux(grid, Point::Ones());
        LevelBoxData<double>    LPhiError(crseLayout, Point::Zeros());
        LevelBoxData<double>    RefluxCorr(crseLayout, Point::Zeros());
        LevelBoxData<double>    RefluxCrse(crseLayout, Point::Zeros());
        LevelBoxData<double>    RefluxFine(crseLayout, Point::Zeros());
        LevelBoxData<double>    LPhiCF(cfLayout, Point::Zeros());
        LevelBoxData<double>    LPhiCFAvg(cfLayout, Point::Zeros());
        LevelBoxData<double>    LPhiSolnCF(cfLayout, Point::Zeros());
        LevelBoxData<double>    LPhiErrorCF(cfLayout, Point::Zeros());
        LevelBoxData<double, DIM> ScaledFlux(crseLayout, Point::Zeros());
        
        auto& PhiCrse = Phi[0]; auto& LPhiCrse = LPhi[0]; auto& FluxCrse = Flux[0];
        auto& PhiFine = Phi[1]; auto& LPhiFine = LPhi[1]; auto& FluxFine = Flux[1];

        Phi.initConvolve(dx[0], f_phi, kx, ky);
        LPhiSoln.initConvolve(dx[0], f_lphi, kx, ky);
        LPhiSoln[0].copyTo(LPhiSolnCF);
        LPhi.setToZero();
        Flux.setToZero();
        LPhiError.setToZero();
        RefluxCorr.setToZero();
        RefluxCrse.setToZero();
        RefluxFine.setToZero();

        h5.writeAMRData(dx[0], Phi, "Phi");

        // BUILD OPERATORS
        std::vector<Stencil<double>> Grad(DIM);
        std::vector<Stencil<double>> Div(DIM);
        for (int dir = 0; dir < DIM; dir++)
        {
            Div[dir]  = 1.0*Shift::Basis(dir, 1) - 1.0*Shift::Zeros();
            Grad[dir] = Stencil<double>::DiffCellToFace(dir);
        }

        // COMPUTE FLUXES AND SOLUTION
        LevelFluxRegister<double> frFine(crseLayout, fineLayout, refRatioVect);
        LevelFluxRegister<double> frCrse(crseLayout, fineLayout, refRatioVect);
        LevelFluxRegister<double> frCorr(crseLayout, fineLayout, refRatioVect);
        
        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            auto& phi_i  = PhiCrse[*iter];
            auto& lphi_i = LPhiCrse[*iter];
            auto& flux_i = FluxCrse[*iter];
            for (int dir = 0; dir < DIM; dir++)
            {
                auto flux_id = slice(flux_i, dir);
                flux_id |= Grad[dir](phi_i, 1.0/dx[0]);
                frCrse.incrementCoarse(flux_id, *iter, 1.0, dir);
                frCorr.incrementCoarse(flux_id, *iter, 1.0, dir);
                lphi_i += Div[dir](flux_id, 1.0/dx[0]);
            }
        }
        
        for (auto iter = fineLayout.begin(); iter.ok(); ++iter)
        {
            auto& phi_i  = PhiFine[*iter];
            auto& lphi_i = LPhiFine[*iter];
            auto& flux_i = FluxFine[*iter];
            for (int dir = 0; dir < DIM; dir++)
            {
                auto flux_id = slice(flux_i, dir);
                flux_id |= Grad[dir](phi_i, 1.0/dx[1]);
                frFine.incrementFine(flux_id, *iter, 1.0, dir);
                frCorr.incrementFine(flux_id, *iter, 1.0, dir);
                lphi_i += Div[dir](flux_id, 1.0/dx[1]);
            }
        }
        LPhiCrse.copyTo(LPhiCF);
        std::cout << "Integral of Phi:              " << Phi.integrate(dx[0])      << std::endl;
        std::cout << "Integral of LPhiSoln (Valid): " << LPhiSoln.integrate(dx[0]) << std::endl;
        double refinedIntegral = LPhiSoln[1].integrate(dx[1]);
        std::cout << "Integral of LPhiSoln (Fine):  " << refinedIntegral << std::endl;
        
        std::cout << std::endl;
        std::cout << "Integral of LPhi Coarse (No Average, No Reflux):    " << LPhi.integrate(dx[0]) << std::endl;
        LPhi.averageDown();
        LPhiCrse.copyTo(LPhiCFAvg);
        std::cout << "Integral of LPhi Coarse (Averaged Down, No Reflux): " << LPhi.integrate(dx[0]) << std::endl;
        
        frFine.reflux(RefluxFine, 1.0/dx[0]);
        frCrse.reflux(RefluxCrse, 1.0/dx[0]);
        frCorr.reflux(RefluxCorr, 1.0/dx[0]);
        frCorr.reflux(LPhiCrse, 1.0/dx[0]);
        
        std::cout << "Integral of LPhi Coarse (Averaged Down, Refluxed):  " << LPhi.integrate(dx[0]) << std::endl;
        std::cout << std::endl;
        std::cout << "Integral of Coarse Reflux / cdx: " << RefluxCrse.integrate(dx[0]) << std::endl;
        std::cout << "\tError: " << RefluxCrse.integrate(dx[0]) - refinedIntegral << std::endl;
        std::cout << "Integral of Fine Reflux / cdx:   " << RefluxFine.integrate(dx[0]) << std::endl;
        std::cout << "\tError: " << -RefluxFine.integrate(dx[0]) - refinedIntegral << std::endl;
        std::cout << std::endl;
        std::cout << "Integral of LPhi Refined Region (Coarse, Before Avg Down): " << LPhiCF.integrate(dx[0]) << std::endl;
        std::cout << "\tError: " << LPhiCF.integrate(dx[0]) - refinedIntegral << std::endl;
        std::cout << "Integral of LPhi Refined Region (Coarse, After  Avg Down): " << LPhiCFAvg.integrate(dx[0]) << std::endl;
        std::cout << "\tError: " << LPhiCFAvg.integrate(dx[0]) - refinedIntegral << std::endl;
        std::cout << "Integral of LPhiSoln Refined Region (Coarse): " << LPhiSolnCF.integrate(dx[0]) << std::endl;
        std::cout << "\tError: " << LPhiSolnCF.integrate(dx[0]) - refinedIntegral << std::endl;
        
        Flux[0].exchange();
        // COMPUTE ERROR IN LPHI AND FLUX
        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            auto& err_i  = LPhiError[*iter];
            auto& lphi_i = LPhi[0][*iter];
            auto& soln_i = LPhiSoln[0][*iter];
            lphi_i.copyTo(err_i);
            err_i -= soln_i;
        }
       
        for (auto iter = cfLayout.begin(); iter.ok(); ++iter)
        {
            auto& err_cf_i  = LPhiErrorCF[*iter];
            auto& lphi_cf_i = LPhiCF[*iter];
            auto& soln_cf_i = LPhiSolnCF[*iter];
            lphi_cf_i.copyTo(err_cf_i);
            err_cf_i -= soln_cf_i;
        }
         
        err[nn] = LPhiError.absMax();    
        err_cf[nn] = LPhiErrorCF.absMax();
        std::cout << std::endl;
        std::cout << "Error in LPhi: " << err[nn] << std::endl;
        std::cout << "Error in LPhiCF: " << err_cf[nn] << std::endl;
        h5.writeAMRData(dx[0], Flux,       "Flux_N%i", nn);
        h5.writeAMRData(dx[0], LPhi,       "LPhi_N%i", nn);
        h5.writeAMRData(dx[0], LPhiSoln,   "LPhiSoln_N%i", nn);
        h5.writeLevel(dx[0],   RefluxCrse, "RefluxCrse_N%i", nn);
        h5.writeLevel(dx[0],   RefluxFine, "RefluxFine_N%i", nn);
        h5.writeLevel(dx[0],   RefluxCorr, "RefluxCorr_N%i", nn);
        h5.writeLevel(dx[0],   LPhiError,  "LPhiError_N%i", nn);
        
        domainSize *= 2;
        std::cout << std::endl;
        std::cout << std::setfill('=') << std::setw(75) << "=" << std::endl << std::endl;; 
    }

    for (int ii = 1; ii < numIter; ii++)
    {
        std::cout << "Convergence Rate: " << log(err[ii-1]/err[ii])/log(2.0) << std::endl;
        std::cout << "Convergence Rate (Refined Region): " << log(err_cf[ii-1]/err_cf[ii])/log(2.0) << std::endl;
    }
}
