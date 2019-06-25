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

double sinX(Proto::Point a_pt, double a_dx)
{
    double x = a_pt[0]*a_dx + a_dx/2.0;
    return sin(x);
}
double cosAvg(Proto::Point a_pt, double a_dx)
{
    double x0 = a_pt[0]*a_dx;
    double x1 = x0 + a_dx;
    return (sin(x1) - sin(x0))/a_dx;
}

double sinAvg(Proto::Point a_pt, double a_dx)
{
    double x0 = a_pt[0]*a_dx;
    double x1 = x0 + a_dx;
    double shift = M_PI/4.0;
    return (-cos(x1-shift) + cos(x0-shift))/a_dx;
}

double sinXcosY(Proto::Point a_pt, double a_dx)
{
    double x0 = a_pt[0]*a_dx;
    double x1 = x0 + a_dx;
    double y0 = a_pt[1]*a_dx;
    double y1 = y0 + a_dx;
    return (-cos(x1) + cos(x0))/a_dx + (sin(y1) - sin(y0))/a_dx;
}

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

#if RUN_TEST==_TRIG_
    std::cout << "Running Simple Trig Test" << std::endl;
#elif RUN_TEST==_DIPOLE_
    std::cout << "Running Dipole Test" << std::endl;
#endif

    double L = 2.0*M_PI;
    
    double DX[NUM_LEVELS];
    for (int ii = 0; ii < NUM_LEVELS; ii++)
    {
        DX[ii] = L/(domainSize*pow(AMR_REFRATIO,ii));
    }

    std::vector<Proto::Box> domainBoxes;
    domainBoxes.resize(NUM_LEVELS);    
    std::vector<Proto::Box> coarseBoxes;
    coarseBoxes.resize(NUM_LEVELS);    
    for (int ii = 0; ii < NUM_LEVELS; ii++)
    {
        int s = ipow(AMR_REFRATIO,ii);
        Proto::Box domain;
#if NESTING==_FULL_NESTING_
        domain = Proto::Box::Cube(domainSize/s);
        if (ii > 0)
        {
            domain = domain.shift(Proto::Point::Ones(domainSize*0.5*(1-1.0/s)));
        }
        domain = domain.refine(s);
#elif NESTING==_X_NESTING_
        domain = Proto::Box(Proto::Point(domainSize/s-1, domainSize-1));
        if (ii > 0)
        {
            domain = domain.shift(Proto::Point::Basis(0, domainSize*0.5*(1-1.0/s)));
        }
        domain = domain.refine(s);
#elif NESTING==_Y_NESTING_
        domain = Proto::Box(Proto::Point(domainSize-1, domainSize/s-1));
        if (ii > 0)
        {
            domain = domain.shift(Proto::Point::Basis(1, domainSize*0.5*(1-1.0/s)));
        }
        domain = domain.refine(s);
#elif NESTING==_NO_NESTING_
        domain = Proto::Box::Cube(domainSize);
        domain = domain.refine(s);
#else
        std::cout << "Could not identify nesting flag: " << NESTING << std::endl;
        return 1;
#endif
        domainBoxes[ii] = domain;
        coarseBoxes[ii] = domain.coarsen(AMR_REFRATIO);
    }
    AMRLayout Layout(domainBoxes, Proto::Point::Ones());
    AMRData<OP::numcomps()> Phi(Layout, OP::ghost(), DX[0], true);
    AMRData<OP::numcomps()> Sln(Layout, OP::ghost(), DX[0], false);
    AMRData<OP::numcomps()> Rhs(Layout, Proto::Point::Zeros(), DX[0], false);
    AMRData<OP::numcomps()> RhsSrc(Layout, Proto::Point::Ones(2), DX[0], false);
    AMRData<OP::numcomps()> RhsTmp(Layout, Proto::Point::Ones(), DX[0], false);
    AMRData<OP::numcomps()> Res(Layout, Proto::Point::Zeros(), DX[0], true);
    AMRData<OP::numcomps()> Err(Layout, Proto::Point::Zeros(), DX[0], true);
    AMRData<OP::numcomps()> Crs(CoarseLayout, Proto::Point::Zeros(), DX[0]*AMR_REFRATIO, false);
#if RUN_TEST==_TRIG_    
    RhsTmp.initialize(
    [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_data, Real a_dx)
    {
        a_data(0) = sinAvg(a_pt, a_dx); 
    }, DX[0]);
    
#elif RUN_TEST==_DIPOLE_
    RhsSrc.initialize(
    [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_data, Real a_dx)
    {
        a_data(0) = sinAvg(a_pt, a_dx); 
        std::vector<double> x0;
        std::vector<double> x1;
        x0.resize(DIM);
        x1.resize(DIM);
        for (int ii = 0; ii < DIM; ii++)
        {
            if (ii == 1) {x0[ii] = 3*L/4.0; x1[ii] = L/4.0;}
            else {
                x0[ii] = L/2.0; x1[ii] = L/2.0;
            }
        }
        a_data(0)  = charge(a_pt, a_dx, x0, L/8.0, 6);
        a_data(0) -= charge(a_pt, a_dx, x1, L/8.0, 6);

    }, DX[0]);

    AMRFAS<LaplaceOp<DATA>> laplaceOp(Layout, DX[NUM_LEVELS-1], 1);
    laplaceOp(RhsTmp, RhsSrc);
    RhsTmp.add(RhsSrc, 1.0/24.0);
#endif    
    
#if RUN_TEST==_TRIG_
    Sln.initialize(
    [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_data, Real a_dx)
    {
        a_data(0) = -sinAvg(a_pt, a_dx); 
    }, DX[0]);
#endif    

RhsTmp.write("RhsTmp.hdf5");
RhsSrc.write("RhsSrc.hdf5");

#if OPERATOR==_MEHRSTELLEN_
    AMRFAS<MehrstellenCorrectionOp<DATA>> correctOp(Layout, DX[NUM_LEVELS-1], 1);
    correctOp(Rhs, RhsTmp);
    Rhs.add(RhsTmp);
#else
    RhsTmp.copyTo(Rhs);
#endif

    cout << "Running full AMR test" << endl;
    AMRFAS<OP> amr_op(Layout, DX[NUM_LEVELS-1], log2(1.0*domainSize) - 1);
    amr_op.residual(Res, Phi, Rhs);
    cout << "Integral of initial Rhs: " << Rhs.integrate() << endl;
    Rhs.write("AMR_Rhs.hdf5");
    for (int nn = 0; nn < numIter; nn++)
    {
        Res.write("AMR_Res_%i.hdf5", nn);
        Phi.write("AMR_Phi_%i.hdf5", nn);
        amr_op.vcycle(Phi, Rhs, Res, nn);
        cout << "Residual: Max = " << scientific << Res.absMax();
        cout << "\t\tIntegral: " << Res.integrate() << endl;
    }
    Res.write("AMR_Res_%i.hdf5", numIter);
    Phi.write("AMR_Phi_%i.hdf5", numIter);
    
    double phiAvg = Phi.integrate() / pow(L,DIM);
    
    std::cout << "Integral of solution: " << Phi.integrate() << std::endl;
    std::cout << "Average value of solution: " << phiAvg << std::endl;
    
#if RUN_TEST==_TRIG_
    Phi.copyTo(Err);
    Err.add(Sln, -1);
    Err.write("AMR_Err.hdf5");
    cout << "Error: " << Err.absMax() << endl;
    
    AMRData<OP::numcomps()> Tmp(Layout, Proto::Point::Zeros(), DX[0], true);
    Tmp.initialize(
            [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_data, Real a_dx)
            {
            a_data(0) = -phiAvg;
            }, DX[0]);
    Phi.add(Tmp);
    std::cout << "Integral of corrected solution: " << Phi.integrate() << std::endl;
    Phi.copyTo(Err);
    Err.add(Sln, -1);
    Err.write("AMR_Err_corrected.hdf5");
    cout << "Corrected Error: " << Err.absMax() << endl;

    AMRData<OP::numcomps()> Eps(Layout, OP::ghost(), DX[0], true);
    AMRData<OP::numcomps()> LEps(Layout, Proto::Point::Zeros(), DX[0], true);
    Phi.copyTo(Eps);
    Eps.add(Sln, -1);
    amr_op(LEps, Eps);
    Eps.write("AMR_Eps.hdf5");
    LEps.write("AMR_Leps.hdf5");
#endif
}
