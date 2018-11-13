#include "LevelData.H"
#include "Multigrid.H"
#include "AMRFAS.H"
#include "Proto.H"
#include "TestOp.H" //definition of OP
#include "UnitTest.H"

using namespace Proto;

int main(int argc, char** argv)
{

    typedef BoxData<Real, NUMCOMPS> BD;
    typedef TestOp<FArrayBox> OP;
    typedef FArrayBox DATA;
    
    if (argc == 2)
    {
      VERBO = atoi(argv[1]);
    } else {
      VERBO = 0;
    }

    cout << "What would you like to run?" << endl;
    cout << "\tTest Set 1: AMR Functions" << endl;
    cout << "\tTest Set 2: Multigrid" << endl;
    cin >> TEST;

    //====================================================================
    BEGIN_TEST_SUITE(1, "AMR Functions");
    BEGIN_TEST("Boundary Condition Interpolation");
    int domainSize = 32;
    Real dx = 2.0*M_PI/domainSize;
    
    Box domainBoxC = Box::Cube(domainSize);
    Box domainBoxF = Box::Cube(domainSize/2).shift(Point::Ones(domainSize/4));
    domainBoxF = domainBoxF.refine(AMR_REFRATIO);
    DisjointBoxLayout coarseLayout, fineLayout;
    buildLayout(coarseLayout, domainBoxC);
    buildLayout(fineLayout, domainBoxF, Point::Zeros());

    LevelData<DATA> LDF(fineLayout, 1, Point::Ones());
    LevelData<DATA> LDC(coarseLayout, 1, Point::Ones());

    auto fiter = fineLayout.dataIterator();
    auto citer = coarseLayout.dataIterator();

    __OUT(2)
    cout << "Coarse Layout" << endl;
    for (citer.begin(); citer.ok(); ++citer)
    {
      cout << Box(coarseLayout[citer]) << endl;
    }

    cout << "Fine Layout" << endl;
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
      cout << Box(fineLayout[fiter]) << endl;
    }
    OUT__

    for (fiter.begin(); fiter.ok(); ++fiter)
    {
      BoxData<double> bd = LDF[fiter];
      bd.setVal(1337);
    }
    for (citer.begin(); citer.ok(); ++citer)
    {
      BoxData<double> bd = LDC[citer];
      bd.setVal(17);
    }
   
    TestOp<DATA> op(fineLayout,dx);
    op.interpBoundary(LDF, LDC);
    
    END_TEST();
    END_TEST_SUITE();

    //====================================================================
    BEGIN_TEST_SUITE(2, "Multigrid");
    BEGIN_TEST("Multigrid");
    int domainSize = 8;
    int numLevels = 5;
    Real dx = 2.0*M_PI/domainSize;
    
    Box domainBox = Box::Cube(domainSize);
    DisjointBoxLayout layout;
    buildLayout(layout, domainBox);

    LevelData<DATA> U(layout, NUMCOMPS, IntVect::Unit);
    LevelData<DATA> F(layout, NUMCOMPS, IntVect::Zero);
    LevelData<DATA> R(layout, NUMCOMPS, IntVect::Zero);
    LevelData<DATA> S(layout, NUMCOMPS, IntVect::Zero);
    
    OP::initialCondition(U,dx);
    OP::forcing(F,dx);
    OP::solution(S,dx);
    
    Multigrid<OP, DATA> mg(layout, dx, numLevels-1);
    int numIter = 20;
    double resnorm = 0.0;
    char fileName[100];
    char fileNameU[100];
    TestOp<FArrayBox> op(layout,dx);
    fileNum = 0;
    for (int ii = 0; ii < numIter; ii++)
    {
        mg.vcycle(U,F); 
        resnorm = op.residual(R,U,F);
        std::cout << "iteration number = " << ii << ", Residual norm: " << resnorm << std::endl;
        sprintf(fileName,"ResV.%i.hdf5",fileNum);
        sprintf(fileNameU,"ResU.%i.hdf5",fileNum);
        writeLevelname(&R,fileName);
        writeLevelname(&U,fileNameU);
        fileNum++;
    }
    auto iter = layout.dataIterator();
    double umax = 0.0;
    double umin = 0.0;
    double error = 0.0;
    for (iter.begin(); iter.ok(); ++iter)
    {
        BoxData<double> s = S[iter()];
        BoxData<double> u = U[iter()];
        umax = max(umax,u.max());
        umin = min(umin,u.min());
        s -= u;
        error = max(s.absMax(),error);
    }
    cout << "Error: " << error << endl;
    cout << "Max: " << umax << " Min: " << umin << endl;
    END_TEST();
    END_TEST_SUITE();

    //==================================================================== 
    BEGIN_TEST_SUITE(3,"Multigrid Operations");
    //==================================================================== 
    /*
    DisjointBoxLayout tempLayout;
    coarsen_dbl(tempLayout, layout,2);

    Box b = Box::Cube(domainSize/2);
    DisjointBoxLayout coarseLayout;
    buildLayout(coarseLayout, b);
    BEGIN_TEST("TestOp::coarsen");
    
    LevelData<DATA> V(layout, NUMCOMPS, IntVect::Unit);
    LevelData<DATA> VC(coarseLayout, NUMCOMPS, IntVect::Unit);
    LevelData<DATA> VC0(coarseLayout, NUMCOMPS, IntVect::Unit);
    
    OP::initialCondition(V,dx);
    TestOp<FArrayBox> op(layout,dx);
    op.coarsen(VC,V);
    VC.copyTo(VC0);
    
    auto citer = coarseLayout.dataIterator();
    auto fiter = layout.dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        for (citer.begin(); citer.ok(); ++citer)
        {
            Box b = layout[fiter()];
            Box bc = coarseLayout[citer()];
            if (bc.refine(MG_REFRATIO).contains(b))
            {
                BoxData<double> v = V[fiter()];
                BoxData<double> vc = VC[citer()];
                for (auto biter = b.begin(); biter != b.end(); ++biter)
                {
                    Point q = *biter;
                    Point p = q/MG_REFRATIO;
                    if ((q % Point::Ones(MG_REFRATIO)) == Point::Zeros())
                    {
                        UNIT_TEST((abs(vc(p) - (v(q) + v(q + Point::Basis(0)) +
                                           v(q+Point::Basis(1)) +
                                           v(q+Point::Ones()))/4.0) < 1e-15));
                    }
                }
            }
        }
    }
    END_TEST();
    //==================================================================== 
    BEGIN_TEST("TestOp::residual");
    
    LevelData<DATA> V(layout, NUMCOMPS, IntVect::Unit);
    LevelData<DATA> F(layout, NUMCOMPS, IntVect::Unit);
    LevelData<DATA> R(layout, NUMCOMPS, IntVect::Unit);
    
    auto iter = layout.dataIterator();
    for (iter.begin(); iter.ok(); ++iter)
    {
        BoxData<double> v = V[iter()];
        BoxData<double> f = F[iter()];
        BoxData<double> r = R[iter()];
        forallInPlace_p([dx](Point p, Var<double>& v, Var<double>& f, Var<double>& r)
        {
            double x = p[0]*dx;
            r(0) = 0.0;
            v(0) = 0.0;
            f(0) = cos(x);
        },v,f,r);
    }
     
    TestOp<FArrayBox> op(layout,dx);
    op.residual(R,V,F);
    
    auto fiter = layout.dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {

    }
    END_TEST();
    //==================================================================== 

    BEGIN_TEST("TestOp::fineCorrection");
   
    LevelData<DATA> VC0(coarseLayout, NUMCOMPS, IntVect::Unit);
    LevelData<DATA> VC(coarseLayout, NUMCOMPS, IntVect::Unit);
    LevelData<DATA> VC_save(coarseLayout, NUMCOMPS, IntVect::Unit);
    LevelData<DATA> V(layout, NUMCOMPS, IntVect::Unit);
    LevelData<DATA> U(layout, NUMCOMPS, IntVect::Unit);
    LevelData<DATA> U_save(layout, NUMCOMPS, IntVect::Unit);

    double dx_c = MG_REFRATIO*dx;
    auto citer = coarseLayout.dataIterator();
    for (citer.begin(); citer.ok(); ++citer)
    {
        Box b = coarseLayout[citer()];
        BoxData<double> vc = VC[citer()];
        BoxData<double> vc0 = VC0[citer()];
        forallInPlace_p([dx_c](Point p, Var<double>& vc, Var<double>& vc0)
        {
            double x = p[0]*dx_c;
            vc0(0) = 0.9*cos(x);
            vc(0) = cos(x);
        },b,vc,vc0);
    }
    VC.copyTo(VC_save);
     
    auto fiter = layout.dataIterator();
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        Box b = layout[fiter()];
        BoxData<double> v = V[fiter()];
        BoxData<double> u = U[fiter()];
        forallInPlace_p([dx](Point p, Var<double>& v, Var<double>& u)
        {
            double x = p[0]*dx;
            u(0) = 0.5*cos(x);
            v(0) = 0.6*cos(x);
        },b,v,u);
    }
    U.copyTo(U_save);

    TestOp<FArrayBox> op(layout,dx);
    op.fineCorrection(U,VC,VC0);
    
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
        for (citer.begin(); citer.ok(); ++citer)
        {
            Box b = layout[fiter()];
            Box bc = coarseLayout[citer()];
            if (bc.refine(MG_REFRATIO).contains(b))
            {
                BoxData<double> u_corr = U[fiter()];
                BoxData<double> u = U_save[fiter()];
                BoxData<double> corr = VC[citer()];
                BoxData<double> vc = VC_save[citer()];
                BoxData<double> vc0 = VC0[citer()];
                __OUT(2) {    
                    cout << "FINE DATA" << endl;
                    u.printData(b);
                    cout << "CORRECTED FINE DATA" << endl;
                    u_corr.printData(b);
                    cout << "CORRECTION" << endl;
                    corr.printData(bc);
                    cout << "UPDATED COARSE DATA" << endl;
                    vc.printData(bc);
                    cout << "SAVED COARSE DATA" << endl;
                    vc0.printData(bc);
                } OUT__
                for (auto biter = b.begin(); biter != b.end(); ++biter)
                {
                    Point q = *biter;
                    Point p = q/MG_REFRATIO;
                    if ((q % Point::Ones(MG_REFRATIO)) == Point::Zeros())
                    {
                        Box K(Point::Ones());
                        for (auto kiter = K.begin(); kiter != K.end(); ++kiter)
                        {
                            __OUT(2) {
                                double error = abs(u_corr(q+(*kiter)) - (vc(p) - vc0(p) + u(q+(*kiter))));
                                if (error > 1e-15)
                                {
                                    cout << "Source -> Dest | " << p << " -> " << q+(*kiter) << " ";
                                    cout << "Error: " << scientific << (u_corr(q+(*kiter)) - (vc(p) - vc0(p) + u(q+(*kiter)))) << " ";
                                    cout << "\t\t Value " << u_corr(q+(*kiter)) << " should be " << (vc(p) - vc0(p) + u(q+(*kiter))) << endl;
                                }
                            } OUT__
                            UNIT_TEST((abs(u_corr(q+(*kiter)) - (vc(p) - vc0(p) + u(q+(*kiter)))) < 1e-15));
                        }
                    }
                }

            }
        }
    }

    END_TEST();
    //==================================================================== 
    */
    END_TEST_SUITE();
}
