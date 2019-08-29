
#include "LevelData.H"
#include "Proto.H"

#define _LAPLACE_ 0
#define _MEHRSTELLEN_ 1

#define _NO_NESTING_ 0
#define _X_NESTING_ 1
#define _Y_NESTING_ 2
#define _FULL_NESTING_ 3

#define _TRIG_ 0
#define _DIPOLE_ 1

#ifndef NESTING
    #define NESTING _X_NESTING_
#endif

#ifndef OPERATOR
    //#define OPERATOR _LAPLACE_
    #define OPERATOR _MEHRSTELLEN_
#endif

#ifndef RUN_TEST
    #define RUN_TEST _TRIG_
#endif

//#define GSRB TRUE

#define NUM_LEVELS 2

#include "AMRFAS.H"
#include "BaseOp.H"
#include "LaplaceOp.H"
#include "MehrstellenOp.H"
#include "MehrstellenCorrectionOp.H"

using namespace std;

double charge(Proto::Point a_pt, double a_dx, std::vector<double> a_x0, double a_r0, int a_p)
{
    std::vector<double> x;
    x.resize(DIM);
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = a_pt[ii]*a_dx + 0.5*a_dx;
    }
    double r = 0.0;
    for (int ii = 0; ii < DIM; ii++)
    {
        r += pow(x[ii] - a_x0[ii],2);
    }
    r = sqrt(r);
    if (r >= a_r0){return 0.0;}
    else {
        return pow(cos(M_PI/2.0*r/a_r0),a_p);
    }
}

int main(int argc, char** argv)
{
    bool output = false;
    int mpi_rank = 0;
#ifdef CH_MPI
    MPI_Init(NULL, NULL);
    int mpi_world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    std::cout << "my rank is: " << mpi_rank << " of " << mpi_world_size << std::endl;
#endif
    char time_table_fname[100];
    sprintf(time_table_fname, "time.table.%i",mpi_rank);
    PR_TIMER_SETFILE(time_table_fname);

#if OPERATOR==_LAPLACE_
    typedef LaplaceOp<FArrayBox> OP;
    if (mpi_rank == 0)
    {
        std::cout << "Using Operator: LaplaceOp" << std::endl;
    }
#elif OPERATOR==_MEHRSTELLEN_
    typedef MehrstellenOp OP;
    if (mpi_rank == 0)
    {
        std::cout << "Using Operator: MehrstellenOp" << std::endl;
    }
#else 
    if (mpi_rank == 0)
    {
        cout << "Could not recognize operator: " << OPERATOR << endl;
    }
    return 1;
#endif
    int numIter = 1;
    int domainSize = 32;
    if (argc >= 2) 
    {
        numIter = atoi(argv[1]);
    }
    if (argc >= 3) {
        domainSize = atoi(argv[2]);
    }

    if (mpi_rank == 0)
    {
        cout << "DomainSize: " << domainSize << endl;
        cout << "Iterations: " << numIter << endl;
    }
#ifdef GSRB
    if (mpi_rank == 0)
    {
        std::cout << "GSRB = TRUE" << std::endl;
    }
#else
    if (mpi_rank == 0)
    {
        std::cout << "GSRB = FALSE" << std::endl;
    }
#endif
    double L = 2.0*M_PI;
    
    int n_runs = 3;
    std::vector<AMRLayout> Layouts;
    Layouts.resize(n_runs);
    if (mpi_rank == 0)
    {
        std::cout << "Defining AMRLayout" << std::endl;
    }
    for (int nn = 0; nn < n_runs; nn++)
    {
        std::vector<Proto::Box> domainBoxes;
        domainBoxes.resize(NUM_LEVELS);    
        for (int ii = 0; ii < NUM_LEVELS; ii++)
        {
            int s = ipow(AMR_REFRATIO,ii);
            Proto::Box domain;
            domain = Proto::Box(Proto::Point(domainSize/s-1, domainSize-1));
            if (ii > 0)
            {
                domain = domain.shift(Proto::Point::Basis(0, domainSize*0.5*(1-1.0/s)));
            }
            domain = domain.refine(s);
            domain = domain.refine(ipow(2,nn));
            domainBoxes[ii] = domain;
            if (mpi_rank == 0) {
                std::cout << domain << ", "; 
            }
        }
        if (mpi_rank == 0)
        {
            std::cout << std::endl;
        }
        Layouts[nn].define(domainBoxes, Proto::Point::Ones());
    }
    if (mpi_rank == 0)
    {
        std::cout << "Building Data Holders" << std::endl;
    }
    std::vector<std::shared_ptr<AMRData<OP::numcomps()>>> Err;
    std::vector<std::shared_ptr<AMRData<OP::numcomps()>>> AllPhi;
    std::vector<std::shared_ptr<AMRData<OP::numcomps()>>> AllRhs;
    Err.resize(n_runs - 1);
    AllPhi.resize(n_runs);
    AllRhs.resize(n_runs);

    if (mpi_rank == 0)
    {
        std::cout << "Executing Convergence Test" << std::endl;
    }
    for (int nn = 0; nn < n_runs; nn++)
    {
        if (mpi_rank == 0)
        {
            std::cout << "\tInitializing Data" << std::endl;
        }
        double DX[NUM_LEVELS];
        for (int ii = 0; ii < NUM_LEVELS; ii++)
        {
            DX[ii] = L/(domainSize*pow(AMR_REFRATIO,ii));
        }
        AMRLayout& Layout = Layouts[nn];
        AllPhi[nn] = std::make_shared<AMRData<OP::numcomps()>>(Layout, OP::ghost(), DX[0], true);
        auto& Phi = *AllPhi[nn];
        AllRhs[nn] = std::make_shared<AMRData<OP::numcomps()>>(Layout, Proto::Point::Zeros(), DX[0], true);
        auto& Rhs = *AllRhs[nn];
        AMRData<OP::numcomps()> RhsSrc(Layout, Proto::Point::Ones(2), DX[0], false);
        AMRData<OP::numcomps()> RhsTmp(Layout, Proto::Point::Ones(), DX[0], false);
        AMRData<OP::numcomps()> Res(Layout, Proto::Point::Zeros(), DX[0], true);
        
        RhsSrc.initialize(
                [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_data, Real a_dx)
                {
                    std::vector<double> x0;
                    std::vector<double> x1;
                    x0.resize(DIM);
                    x1.resize(DIM);
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        if (ii == 1) {x0[ii] = 3*L/4.0; x1[ii] = L/4.0;}
                        else {
                            x0[ii] = L/4.0; x1[ii] = 3.0*L/4.0;
                        }
                    }
                    a_data(0)  = charge(a_pt, a_dx, x0, L/4.0, 6);
                    a_data(0) -= charge(a_pt, a_dx, x1, L/4.0, 6);

                }, DX[0]);

       RhsSrc.toCellAverage(RhsTmp);
        
#if OPERATOR==_MEHRSTELLEN_
        AMRFAS<MehrstellenCorrectionOp> correctOp(Layout, DX[NUM_LEVELS-1], 1);
        correctOp(Rhs, RhsTmp);
        Rhs.increment(RhsTmp);
#else
        RhsTmp.copyTo(Rhs);
#endif
        if (mpi_rank == 0)
        {
            std::cout << "\tBuilding Operator" << std::endl;
        }
    
        AMRFAS<OP> amr_op(Layout, DX[NUM_LEVELS-1], log2(1.0*domainSize) - 1);
        amr_op.residual(Res, Phi, Rhs);

        if (mpi_rank == 0)
        {
            std::cout << "\tComputing Initial Integrals..." << std::endl;
        }
        double rhsInt = Rhs.integrate();
        double resInt = Res.integrate();
        if (mpi_rank == 0)
        {
            cout << "Integral of initial Rhs: " << rhsInt << endl;
            cout << "Integral of initial Res: " << resInt << endl;
        }
        //Rhs.write("AMR_Rhs_N%i.hdf5", nn);
        //Res.write("AMR_Res_N%i_0.hdf5", nn);
        //Phi.write("AMR_Phi_N%i_0.hdf5", nn);
        for (int jj = 0; jj < numIter; jj++)
        {
#ifdef CH_MPI
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            amr_op.vcycle(Phi, Rhs, Res, nn);
            Real resMax = Res.absMax();
            Real resInt = Res.integrate();
            if (mpi_rank == 0)
            {
                cout << "Residual: Max = " << scientific << resMax;
                cout << "\t\tIntegral: " << resInt << endl;
            }
            //Res.write("AMR_Res_N%i_%i.hdf5", nn, jj+1);
            //Phi.write("AMR_Phi_N%i_%i.hdf5", nn, jj+1);
        }
        //Res.write("AMR_Res_N%i.hdf5", nn);
        //Phi.write("AMR_Phi_N%i.hdf5", nn);

        Real phiInt = Phi.integrate();
        double phiAvg = phiInt / pow(L,DIM);
        if (mpi_rank == 0)
        {
            std::cout << "Integral of solution: " << phiInt << std::endl;
            std::cout << "Average value of solution: " << phiAvg << std::endl;
        }

        if (nn > 0) {
        
            Err[nn-1] = make_shared<AMRData<OP::numcomps()>>(Layouts[nn-1], Proto::Point::Zeros(), DX[0]*AMR_REFRATIO, true);
            Phi.coarsenTo(*Err[nn-1]);
            (*Err[nn-1]).increment(*AllPhi[nn-1], -1);
            (*Err[nn-1]).write("Error_%i.hdf5", nn-1);
        }
        domainSize *= 2;
    } //end runs
    for (int ii = 1; ii < n_runs - 1; ii++)
    {
        Real err0 = (*Err[ii-1]).absMax();
        Real err1 = (*Err[ii]).absMax();
        if (mpi_rank == 0)
        {
            std::cout << "Error 1: " << err0 << ", Error 2: " << err1 << std::endl;
            std::cout << "Rate: " << log2(err0/err1) << std::endl;
        }
    }

    PR_TIMER_REPORT();
#ifdef CH_MPI
    MPI_Finalize();
#endif
}
