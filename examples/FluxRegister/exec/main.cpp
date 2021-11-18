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

void f_error(Point& a_pt, Var<double, DIM>& a_error, Var<double>& a_reflux, Var<double, DIM>& a_flux, Box a_cfBox)
{
    for (int dir = 0; dir < DIM; dir++)
    {
        Box upper = a_cfBox.adjacent(dir, Side::Hi, 1);
        Box lower = a_cfBox.adjacent(dir, Side::Lo, 1);
        if (upper.contains(a_pt))
        {
            a_error(dir) = a_reflux(0) - a_flux(dir);
        } else if (lower.contains(a_pt))
        {
            a_error(dir) = a_reflux(0) + a_flux(dir);
        } else {
            a_error(dir) = 0;
        }
    }
}

int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init (&argc, &argv);
    setPoutBaseName("FluxRegister");
#endif

    // READ INPUT PARAMETERS
    HDF5Handler h5;
    
    InputArgs args;
    args.parse(); //assumes a file "inputs" exists.
    int domainSize = args.get("domainSize");
    int boxSize = args.get("boxSize");
    int numIter = args.get("numIter");
    //int refRatio = args.get("refRatio");
    int refRatio = PR_AMR_REFRATIO;
    PR_TIMER_SETFILE(to_string(domainSize) + ".proto.time.table");
    PR_TIMERS("main");

    // NON INPUT PARAMETERS
    //double kx = sqrt(2.0)/2.0;
    //double ky = sqrt(2.0)/2.0;
    double kx = 1;
    double ky = 0;
    Point refRatioVect = Point::Ones(refRatio);
    Point boxSizeVect = Point::Ones(boxSize);
    
    double err[numIter];
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

        // BUILD DATA HOLDERS
        AMRData<double>         Phi(grid, Point::Ones(3));
        AMRData<double>         LPhi(grid, Point::Zeros());
        AMRData<double>         LPhiSoln(grid, Point::Zeros());
        AMRData<double, DIM>    Flux(grid, Point::Ones());
        LevelBoxData<double>    LPhiError(crseLayout, Point::Zeros());
        LevelBoxData<double>    refluxCorr(crseLayout, Point::Zeros());
        LevelBoxData<double>    refluxCrse(crseLayout, Point::Zeros());
        LevelBoxData<double>    refluxFine(crseLayout, Point::Zeros());
        auto& PhiCrse = Phi[0]; auto& LPhiCrse = LPhi[0]; auto& FluxCrse = Flux[0];
        auto& PhiFine = Phi[1]; auto& LPhiFine = LPhi[1]; auto& FluxFine = Flux[1];

        Phi.initConvolve(dx[0], f_phi, kx, ky);
        LPhiSoln.initConvolve(dx[0], f_lphi, kx, ky);
        LPhi.setToZero();
        Flux.setToZero();
        LPhiError.setToZero();
        refluxCorr.setToZero();
        refluxCrse.setToZero();
        refluxFine.setToZero();

        h5.writeAMRData(dx[0], Phi, "Phi");

        // BUILD OPERATORS
        std::vector<Stencil<double>> Grad(DIM);
        std::vector<Stencil<double>> Div(DIM);
        for (int dir = 0; dir < DIM; dir++)
        {
            Div[dir] =  1.0*Shift::Basis(dir, 1) - 1.0*Shift::Zeros();
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
        std::cout << "Integral of Phi: "      << Phi.integrate(dx[0])      << std::endl;
        std::cout << "Integral of LPhiSoln: " << LPhiSoln.integrate(dx[0]) << std::endl;
        std::cout << std::endl;
        std::cout << "Integral of LPhiCrse (No Average, No Reflux): " << LPhi.integrate(dx[0]) << std::endl;
        LPhi.averageDown();
        std::cout << "Integral of LPhiCrse (Averaged Down, No Reflux): " << LPhi.integrate(dx[0]) << std::endl;
        
        frFine.reflux(refluxFine, 1.0);
        frCrse.reflux(refluxCrse, 1.0);
        frCorr.reflux(refluxCorr, 1.0);
        frCorr.reflux(LPhiCrse, 1.0);
        std::cout << "Integral of LPhiCrse (Averaged Down, Refluxed): " << LPhi.integrate(dx[0]) << std::endl;
        
        // COMPUTE ERROR IN LPHI
        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            auto& err_i = LPhiError[*iter];
            auto& lphi_i = LPhi[0][*iter];
            auto& soln_i = LPhiSoln[0][*iter];
            lphi_i.copyTo(err_i);
            err_i -= soln_i;
        }
        
        err[nn] = LPhiError.absMax();    
        std::cout << "Error in LPhi: " << err[nn] << std::endl;

        h5.writeAMRData(dx[0], Flux, "Flux");
        h5.writeAMRData(dx[0], LPhi, "LPhi");
        h5.writeAMRData(dx[0], LPhiSoln, "LPhiSoln");
        h5.writeLevel(dx[0],   refluxCrse, "RefluxCrse");
        h5.writeLevel(dx[0],   refluxFine, "RefluxFine");
        h5.writeLevel(dx[0],   refluxCorr, "RefluxCorr");
        h5.writeLevel(dx[0],   LPhiError,  "LPhiError");
       
        
        /*
        std::cout << std::endl;
        std::cout << "Local Sums: Reflux Correction" << std::endl;
        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            std::cout << "\tSum in box " << iter.box() << ": " << refluxCorr[*iter].sum() << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Local Sums: Coarse Fluxes" << std::endl;
        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            std::cout << "\tSum in box " << iter.box() << ": " << refluxCrse[*iter].sum() << std::endl;
        //    refluxCrse[*iter].printData();
        }
        std::cout << std::endl;
        std::cout << "Local Sums: Fine Fluxes" << std::endl;
        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            std::cout << "\tSum in box " << iter.box() << ": " << refluxFine[*iter].sum() << std::endl;
        }
        */
        domainSize *= 2;
        std::cout << std::endl;
        std::cout << std::setfill('=') << std::setw(75) << "=" << std::endl << std::endl;; 
    }

    for (int ii = 1; ii < numIter; ii++)
    {
        std::cout << "Convergence Rate: " << log(err[ii-1]/err[ii])/log(2.0) << std::endl;
    }
}
