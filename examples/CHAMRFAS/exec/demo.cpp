#include "LevelData.H"
#include "Proto.H"

#define _LAPLACE_ 0
#define _MEHRSTELLEN_ 1

#define _NO_NESTING_ 0
#define _X_NESTING_ 1
#define _Y_NESTING_ 2
#define _FULL_NESTING_ 3

#define NESTING _FULL_NESTING_

//#define OPERATOR _LAPLACE_
#define OPERATOR _MEHRSTELLEN_

#define GSRB TRUE

#define NUM_LEVELS 2

#include "AMRFAS.H"
#include "BaseOp.H"
#include "LaplaceOp.H"
#include "MehrstellenOp.H"
#include "MehrstellenCorrectionOp.H"

using namespace std;

double cosAvg(Proto::Point a_pt, double a_dx)
{
    double x0 = a_pt[0]*a_dx;
    double x1 = x0 + a_dx;
    return (sin(x1) - sin(x0))/a_dx;
}

int main(int argc, char** argv)
{
    typedef FArrayBox DATA;

#if OPERATOR==_LAPLACE_
    typedef LaplaceOp<FArrayBox> OP;
    std::cout << "Using Operator: LaplaceOp" << std::endl;
#elif OPERATOR==_MEHRSTELLEN_
    typedef MehrstellenOp<FArrayBox> OP;
    std::cout << "Using Operator: MehrstellenOp" << std::endl;
#else 
    cout << "Could not recognize operator: " << OPERATOR << endl;
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

    cout << "DomainSize: " << domainSize << endl;
    cout << "Iterations: " << numIter << endl;

#ifdef GSRB
    std::cout << "GSRB = TRUE" << std::endl;
#else
    std::cout << "GSRB = FALSE" << std::endl;
#endif

    double L = 2.0*M_PI;
    double DX[NUM_LEVELS];
    for (int ii = 0; ii < NUM_LEVELS; ii++)
    {
        DX[ii] = L/(domainSize*pow(AMR_REFRATIO,ii));
    }
    
    std::vector<DisjointBoxLayout> Layouts;
    std::vector<DisjointBoxLayout> TempLayouts;

    std::vector<std::shared_ptr<LevelData<DATA>>> Phi;
    std::vector<std::shared_ptr<LevelData<DATA>>> Sln;
    std::vector<std::shared_ptr<LevelData<DATA>>> Rhs;
    std::vector<std::shared_ptr<LevelData<DATA>>> Res;
    std::vector<std::shared_ptr<LevelData<DATA>>> Err;

    Layouts.resize(NUM_LEVELS);
    TempLayouts.resize(NUM_LEVELS-1);
    Phi.resize(NUM_LEVELS);
    Sln.resize(NUM_LEVELS);
    Rhs.resize(NUM_LEVELS);
    Res.resize(NUM_LEVELS);
    Err.resize(NUM_LEVELS);

    for (int ii = 0; ii < NUM_LEVELS; ii++)
    {
        auto& layout = Layouts[ii];
        int s = ipow(AMR_REFRATIO,ii);
#if NESTING==_FULL_NESTING_
        Box domain = Proto::Box::Cube(domainSize/s);
        if (ii > 0)
        {
            domain = domain.shift(Proto::Point::Ones(domainSize*0.5*(1-1.0/s)));
        }
        domain = domain.refine(s);
        if (ii == 0)
        {
            buildLayout(layout, domain, Proto::Point::Ones());
        } else {
            buildLayout(layout, domain, Proto::Point::Zeros());
        }
#elif NESTING==_X_NESTING_
        Box domain = Proto::Box(Proto::Point(domainSize/s-1, domainSize-1));
        if (ii > 0)
        {
            domain = domain.shift(Proto::Point::Basis(0, domainSize*0.5*(1-1.0/s)));
        }
        domain = domain.refine(s);
        if (ii == 0)
        {
            buildLayout(layout, domain, Proto::Point::Ones());
        } else {
            buildLayout(layout, domain, Proto::Point::Basis(1));
        }
#elif NESTING==_Y_NESTING_
        Box domain = Proto::Box(Proto::Point(domainSize-1, domainSize/s-1));
        if (ii > 0)
        {
            domain = domain.shift(Proto::Point::Basis(1, domainSize*0.5*(1-1.0/s)));
        }
        domain = domain.refine(s);
        if (ii == 0)
        {
            buildLayout(layout, domain, Proto::Point::Ones());
        } else {
            buildLayout(layout, domain, Proto::Point::Basis(0));
        }
#elif NESTING==_NO_NESTING_
        Box domain = Proto::Box::Cube(domainSize);
        domain = domain.refine(s);
        buildLayout(layout, domain, Proto::Point::Ones());
#else
        std::cout << "Could not identify nesting flag: " << NESTING << std::endl;
        return 1;
#endif
        if (ii > 0)
        {
            auto tempLayout = TempLayouts[ii-1];
            coarsen_dbl(tempLayout, layout, AMR_REFRATIO);
        }

        Phi[ii] = std::make_shared<LevelData<DATA>>(layout, 1, OP::ghost());
        Sln[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
        Err[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
        Res[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
        Rhs[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
    
        auto& phiLevel = *(Phi[ii]);
        auto& slnLevel = *(Sln[ii]);
        auto& rhsLevel = *(Rhs[ii]);
        auto& resLevel = *(Res[ii]);
        auto& errLevel = *(Err[ii]);

        auto iter = layout.dataIterator();
        for (iter.begin(); iter.ok(); ++iter)
        {
            OP::patch phi = phiLevel[iter];
            OP::patch rhs = rhsLevel[iter];
            OP::patch res = resLevel[iter];
            OP::patch sln = slnLevel[iter];
            OP::patch err = errLevel[iter];

            phi.setVal(0);
            err.setVal(0);
            res.setVal(0);
            double dx = DX[ii]; 
            forallInPlace_p(
                [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_rhs, OP::var& a_sln)
                {
                    a_rhs(0) = cosAvg(a_pt, dx); 
                    a_sln(0) = -cosAvg(a_pt, dx);
                }, rhs, sln);
            forallInPlace_p(
                [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_phi)
                {
                    a_phi(0) = -cosAvg(a_pt, dx); 
                }, phi);
        }
        phiLevel.exchange();
    }

#if OPERATOR==_MEHRSTELLEN_
    std::vector<std::shared_ptr<LevelData<DATA>>> RhsTemp;
    RhsTemp.resize(NUM_LEVELS);
    for (int ii = 0; ii < NUM_LEVELS; ii++)
    {
        RhsTemp[ii] = std::make_shared<LevelData<DATA>>(Layouts[ii], 1, Proto::Point::Ones());
        auto& tmpLevel = *(RhsTemp[ii]);
        auto iter = Layouts[ii].dataIterator();
        double dx = DX[ii];
        for (iter.begin(); iter.ok(); ++iter)
        {
            OP::patch tmp = tmpLevel[iter];
            forallInPlace_p(
                [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_tmp)
                {
                    a_tmp(0) = cosAvg(a_pt, dx);
                }, tmp);
        }
        tmpLevel.exchange();
    }

    AMRFAS<MehrstellenCorrectionOp<DATA>> correctOp(Layouts, DX[NUM_LEVELS-1], 1);
    correctOp(Rhs, RhsTemp);
    
    for (int ii = 0; ii < NUM_LEVELS; ii++)
    {
        auto& tmpLevel = *(RhsTemp[ii]);
        auto& rhsLevel = *(Rhs[ii]);
        auto iter = Layouts[ii].dataIterator();
        for (iter.begin(); iter.ok(); ++iter)
        {
            OP::patch rhs = rhsLevel[iter];
            OP::patch tmp = tmpLevel[iter];
            rhs += tmp;
        }
    }
#endif

    cout << "Running full AMR test" << endl;
    AMRFAS<OP> amr_op(Layouts, DX[NUM_LEVELS-1], log2(1.0*domainSize) - 1);
    amr_op.residual(Res, Phi, Rhs);

    for (int nn = 0; nn < numIter; nn++)
    {
        amr_op.write(Res, "AMR_Res_%i.hdf5", nn);
        amr_op.write(Phi, "AMR_Phi_%i.hdf5", nn);
        amr_op.vcycle(Phi, Rhs, Res, nn);
        cout << "Residual: Max = " << scientific << absMax(Res);
        cout << "\t\tIntegral: " << integrate(Res, DX[0]) << endl;
    }
    amr_op.write(Res, "AMR_Res_%i.hdf5", numIter);
    amr_op.write(Phi, "AMR_Phi_%i.hdf5", numIter);

    for (int ii = 0; ii < NUM_LEVELS; ii++)
    {
        auto& phiLevel = *(Phi[ii]);
        auto& slnLevel = *(Sln[ii]);
        auto& errLevel = *(Err[ii]);
        auto iter = Layouts[ii].dataIterator();
        for (iter.begin(); iter.ok(); ++iter)
        {
            OP::patch sln = slnLevel[iter];
            OP::patch phi = phiLevel[iter];
            OP::patch err = errLevel[iter];
            phi.copyTo(err);
            err -= sln;
        }
    }
    amr_op.write(Err, "AMR_Err.hdf5");
    Real error = absMax(Err, 0, true, 1);

    cout << "Error: " << error << endl;
}
