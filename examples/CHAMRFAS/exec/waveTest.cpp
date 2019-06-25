#include "LevelData.H"
#include "Proto.H"

#define _LAPLACE_ 0
#define _MEHRSTELLEN_ 1

#define _NO_NESTING_ 0
#define _X_NESTING_ 1
#define _Y_NESTING_ 2
#define _FULL_NESTING_ 3

#ifndef NESTING
    #define NESTING _X_NESTING_
#endif

#ifndef OPERATOR
    //#define OPERATOR _LAPLACE_
    #define OPERATOR _MEHRSTELLEN_
#endif

//#define GSRB TRUE

#define NUM_LEVELS 2

#include "AMRFAS.H"
#include "BaseOp.H"
#include "LaplaceOp.H"
#include "MehrstellenOp.H"
#include "MehrstellenCorrectionOp.H"

#define FREQ 1
#define ALPHA 2*FREQ*M_PI

using namespace std;

double wave(Proto::Point a_pt, double a_dx, std::vector<double> a_x0)
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
    if (r >= 1){return 0.0;}
    else {
        return (r - r*r)*pow(sin(FREQ*M_PI*r), 2);
    }
}

double waveSoln(Proto::Point a_pt, double a_dx, std::vector<double> a_x0)
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
    if (r >= 1)
    {
        return (-1.0/210 - 12.0/pow(ALPHA, 4) + 360.0/pow(ALPHA, 6))/r;
    } else {
        return  pow(r,6)/84.0 - pow(r,5)/30.0 + pow(r,4)/40.0 +
                60.0/pow(ALPHA,6) - 9.0/pow(ALPHA,4) -1.0/120 + 120.0/(pow(ALPHA,6)*r) +
                (
                    -120.0/(r*pow(ALPHA, 6)) - 9.0/pow(ALPHA,4) +300.0/pow(ALPHA, 6) + 36.0*r/pow(ALPHA, 4) + pow(r,2)/(2.0*pow(ALPHA, 2)) -
                    30.0*pow(r, 2)/pow(ALPHA, 4) - pow(r, 3)/pow(ALPHA, 2) + pow(r, 4)/(2.0*pow(ALPHA, 2))
                )*cos(ALPHA*r) +
                (
                    12.0/(r*pow(ALPHA, 5)) - 360.0/(r*pow(ALPHA, 7)) - 96.0/pow(ALPHA, 5) + 120.0*r/pow(ALPHA, 5) -
                    3.0*r/pow(ALPHA, 3) + 8.0*pow(r, 2)/pow(ALPHA, 3) - 5.0*pow(r/ALPHA, 3)
                )*sin(ALPHA*r);
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

    double L = 2.0*M_PI;
    
    double DX[NUM_LEVELS];
    for (int ii = 0; ii < NUM_LEVELS; ii++)
    {
        DX[ii] = L/(domainSize*pow(AMR_REFRATIO,ii));
    }

    std::vector<Proto::Box> domainBoxes;
    domainBoxes.resize(NUM_LEVELS);    
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
    }
    AMRLayout Layout(domainBoxes, Proto::Point::Ones());
    AMRData<OP::numcomps()> Phi(Layout, OP::ghost(), DX[0], true);
    AMRData<OP::numcomps()> Sln(Layout, OP::ghost(), DX[0], false);
    AMRData<OP::numcomps()> Rhs(Layout, Proto::Point::Zeros(), DX[0], false);
    AMRData<OP::numcomps()> RhsSrc(Layout, Proto::Point::Ones(2), DX[0], false);
    AMRData<OP::numcomps()> RhsTmp(Layout, Proto::Point::Ones(), DX[0], false);
    AMRData<OP::numcomps()> Res(Layout, Proto::Point::Zeros(), DX[0], true);
    AMRData<OP::numcomps()> Err(Layout, Proto::Point::Zeros(), DX[0], true);

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
                x0[ii] = L/2.0; x1[ii] = L/2.0;
            }
        }
        a_data(0)  = wave(a_pt, a_dx, x0);
        a_data(0) -= wave(a_pt, a_dx, x1);

    }, DX[0]);

    AMRFAS<LaplaceOp<DATA>> laplaceOp(Layout, DX[NUM_LEVELS-1], 1);
    laplaceOp(RhsTmp, RhsSrc);
    RhsTmp.add(RhsSrc, 1.0/24.0);
   
    Sln.initialize(
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
                x0[ii] = L/2.0; x1[ii] = L/2.0;
            }
        }
        a_data(0)  = waveSoln(a_pt, a_dx, x0);
        a_data(0) -= waveSoln(a_pt, a_dx, x1);
        
    }, DX[0]);

RhsTmp.write("RhsTmp.hdf5");
RhsSrc.write("RhsSrc.hdf5");
Sln.write("Soln.hdf5");

#if OPERATOR==_MEHRSTELLEN_
    AMRFAS<MehrstellenCorrectionOp<DATA>> correctOp(Layout, DX[NUM_LEVELS-1], 1);
    correctOp(Rhs, RhsTmp);
    Rhs.add(RhsTmp);
#else
    RhsTmp.copyTo(Rhs);
#endif

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
}
