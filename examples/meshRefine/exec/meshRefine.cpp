#include "Proto.H"

using namespace Proto;

void foo (Point& a_pt, Var<double>& a_data, double a_dx)
{
    double x = a_pt[0]*a_dx + a_dx/2.0;
    double y = a_pt[1]*a_dx + a_dx/2.0;
    
    a_data(0) = sin(2.0*M_PI*(x + y));
}

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    int domainSize = 32;
    int boxSize = 8;
    if (procID() == 0)
    {
        if (argc > 1)
        {
            domainSize = atoi(argv[1]);
        }
        if (argc > 2)
        {
            boxSize = atoi(argv[2]);
        }
    }

#ifdef PR_MPI
    MPI_Bcast(&domainSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&boxSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    double L = 1.0;
    double k = 2.0*M_PI;
    int refRatio = 2;

    for (int nn = 0; nn < 3; nn++)
    {
        double dx = L / domainSize;
        double fdx = dx/2.0;

        int fineDomainSize = domainSize * refRatio;
        Box crseDomain = Box::Cube(domainSize);
        Box fineDomain = Box::Cube(fineDomainSize);
        Box subDomain  = Box::Cube(fineDomainSize / 2);
        //subDomain = subDomain.shift(Point::Ones(fineDomainSize / 4));
        Box subDomainPatches = subDomain.coarsen(boxSize);
        std::vector<Point> subDomainPatchVector;
        for (auto biter = subDomainPatches.begin();
                biter != subDomainPatches.end(); ++biter)
        {
            subDomainPatchVector.push_back(*biter);
            subDomainPatchVector.push_back(*biter + subDomainPatches.high() + Point::Ones());
        }
        
        Point boxSizeVect = Point::Ones(boxSize);
        Point ghost = Point::Ones();
        std::array<bool, DIM> periodicity;
        for (int ii = 0; ii < DIM; ii++) { periodicity[ii] = true; }

        ProblemDomain crseProblemDomain(crseDomain, periodicity);
        ProblemDomain fineProblemDomain(fineDomain, periodicity);

        DisjointBoxLayout crseLayout(crseProblemDomain, boxSizeVect);
        DisjointBoxLayout fineLayout(fineProblemDomain, subDomainPatchVector, boxSizeVect);
        DisjointBoxLayout crseFineLayout = fineLayout.coarsen(Point::Ones(refRatio));

        crseLayout.print("Coarse Layout");
        fineLayout.print("Fine Layout"); 
        crseFineLayout.print("Coarse Fine Layout");

        auto laplace = Stencil<double>::Laplacian();

        LevelBoxData<double> phi_c(crseLayout,  ghost);
        LevelBoxData<double> phi_f(fineLayout,  ghost);
        LevelBoxData<double> Lphi_c(crseLayout, ghost);
        LevelBoxData<double> Lphi_f(fineLayout, ghost);
        LevelBoxData<double> err_c(crseLayout,  Point::Zeros());
        LevelBoxData<double> err_f(fineLayout,  Point::Zeros());

        phi_c.initialize(foo, dx);
        phi_f.initialize(foo, fdx);

        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            auto& phi = phi_c[*iter];
            auto& Lphi = Lphi_c[*iter];

            Lphi |= laplace(phi, -1.0/(2.0*k*k*dx*dx));
        }

        for (auto iter = fineLayout.begin(); iter.ok(); ++iter)
        {
            auto& phi = phi_f[*iter];
            auto& Lphi = Lphi_f[*iter];

            Lphi |= laplace(phi, -1.0/(2.0*k*k*fdx*fdx));
        }

        HDF5Handler h5;
        std::vector<std::string> vars_phi  = {"phi"};
        std::vector<std::string> vars_err  = {"error"};
        std::vector<std::string> vars_Lphi = {"Lphi"};

        h5.writeLevel(Lphi_c, vars_Lphi, dx,   "Lphi_C_N%i", nn);
        h5.writeLevel(Lphi_f, vars_Lphi, dx/2, "Lphi_F_N%i", nn);
        h5.writeLevel(phi_c,  vars_phi,  dx,   "phi_C_N%i", nn);
        h5.writeLevel(phi_f,  vars_phi,  dx/2, "phi_F_N%i", nn);

        for (auto iter = crseLayout.begin(); iter.ok(); ++iter)
        {
            auto& err = err_c[*iter];
            auto& phi = phi_c[*iter];
            auto& Lphi = Lphi_c[*iter];
            Lphi.copyTo(err);
            err -= phi;
        }

        for (auto iter = fineLayout.begin(); iter.ok(); ++iter)
        {
            auto& err = err_f[*iter];
            auto& phi = phi_f[*iter];
            auto& Lphi = Lphi_f[*iter];
            Lphi.copyTo(err);
            err -= phi;
        }

        h5.writeLevel(err_c, vars_err, dx, "Error_C_N%i", nn);
        h5.writeLevel(err_f, vars_err, dx/2, "Error_F_N%i", nn);
    
        domainSize *= 2;
    }
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
