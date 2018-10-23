#include "LevelData.H"
#include "Multigrid.H"
#include "AMRFAS.H"
#include "Proto.H"
#include "TestOp.H" //definition of OP
#include "NeighborIterator.H"
#include "CFRegion.H"

extern int fileNum;
#ifdef UNIT_TESTING
#include "UnitTest.H"
#endif

using namespace Proto;

int main(int argc, char** argv)
{

    #ifdef UNIT_TESTING
    _TEST_ = 1;
    if (argc == 2)
    {
        _TEST_ = atoi(argv[1]);
        _VERBO_ = 1; 
    }
    if (argc == 3)
    {
        _TEST_ = atoi(argv[1]);
        _VERBO_ = atoi(argv[2]); 
    }
    #endif

    typedef BoxData<Real, NUMCOMPS> BD;
    typedef TestOp<FArrayBox> OP;
    typedef FArrayBox DATA;

    int domainSize = 32;
    int numLevels = 7;
    Real dx = 2.0*M_PI/domainSize;
   
    bool periodic[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        periodic[ii] = true;
    }
    
    Bx domainBoxC = Bx::Cube(domainSize);
    Bx domainBoxF = Bx::Cube(domainSize/2);
    domainBoxF = domainBoxF.refine(AMR_REFRATIO);
    DisjointBoxLayout coarseLayout, fineLayout;
    buildLayout(coarseLayout, domainBoxC);
    buildLayout(fineLayout, domainBoxF, Point::Zeros());

    auto citer = coarseLayout.dataIterator();
    auto fiter = fineLayout.dataIterator();

    cout << "Coarse Layout: " << endl;
    for (citer.begin(); citer.ok(); ++citer)
    {
      cout << Bx(coarseLayout[citer]) << endl;
    }

    cout << "Fine Layout: " << endl;
    for (fiter.begin(); fiter.ok(); ++fiter)
    {
      cout << Bx(fineLayout[fiter]) << endl;
    }

    int patchSizeC = MAXBOXSIZE;
    int patchSizeF = MAXBOXSIZE;
    for (int ii = 0; ii < 0; ii++)
    {
      patchSizeC = min((int)domainBoxC.size(ii),patchSizeC);
      patchSizeF = min((int)domainBoxF.size(ii),patchSizeF);
    }
    Bx bitboxC = domainBoxC.coarsen(patchSizeC);
    Bx bitboxF = domainBoxF.coarsen(patchSizeF);
    Bx bitboxFC = bitboxF.coarsen(AMR_REFRATIO);
    Bx boundsBox = bitboxFC.grow(1);
    cout << "Coarse Bitbox: " << bitboxC << endl;
    cout << "Fine Bitbox: " << bitboxF << endl;
    cout << "Coarsened Fine Bitbox: " << bitboxFC << endl;
    std::set<Point> bdryF;
    std::set<Point> bdryC;
    
    for (int ii = 0; ii < DIM; ii++)
    {
      Bx lo = (boundsBox.flatten(ii,false) & bitboxC);
      for (auto iter = lo.begin(); iter != lo.end(); ++iter)
      {
        bdryC.insert(*iter);
      }
      Bx hi = (boundsBox.flatten(ii,true) & bitboxC);
      for (auto iter = hi.begin(); iter != hi.end(); ++iter)
      {
        bdryC.insert(*iter);
      }
    }

    cout << "Boundary: " << endl;
    for (auto iter = bdryC.begin(); iter != bdryC.end(); ++iter)
    {
      Bx b(*iter,*iter);
      b = b.refine(min((int)domainBoxC.size(0),MAXBOXSIZE));
      cout << b << endl;
    }
    cout << endl;

    #ifndef UNIT_TESTING
   
    /* Multigrid test
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

    */
    #else
    DisjointBoxLayout tempLayout;
    coarsen_dbl(tempLayout, layout,2);

    Box b = Bx::Cube(domainSize/2);
    DisjointBoxLayout coarseLayout;
    buildLayout(coarseLayout, b);
    //==================================================================== 
    BEGIN_TEST_SUITE(1,"TestOp");
    //==================================================================== 
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
            Bx b = layout[fiter()];
            Bx bc = coarseLayout[citer()];
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
        Bx b = coarseLayout[citer()];
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
        Bx b = layout[fiter()];
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
            Bx b = layout[fiter()];
            Bx bc = coarseLayout[citer()];
            if (bc.refine(MG_REFRATIO).contains(b))
            {
                BoxData<double> u_corr = U[fiter()];
                BoxData<double> u = U_save[fiter()];
                BoxData<double> corr = VC[citer()];
                BoxData<double> vc = VC_save[citer()];
                BoxData<double> vc0 = VC0[citer()];
                if (_VERBO_ > 1)
                {
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
                }
                for (auto biter = b.begin(); biter != b.end(); ++biter)
                {
                    Point q = *biter;
                    Point p = q/MG_REFRATIO;
                    if ((q % Point::Ones(MG_REFRATIO)) == Point::Zeros())
                    {
                        Bx K(Point::Ones());
                        for (auto kiter = K.begin(); kiter != K.end(); ++kiter)
                        {
                            if (_VERBO_ > 1)
                            {
                                double error = abs(u_corr(q+(*kiter)) - (vc(p) - vc0(p) + u(q+(*kiter))));
                                if (error > 1e-15)
                                {
                                    cout << "Source -> Dest | " << p << " -> " << q+(*kiter) << " ";
                                    cout << "Error: " << scientific << (u_corr(q+(*kiter)) - (vc(p) - vc0(p) + u(q+(*kiter)))) << " ";
                                    cout << "\t\t Value " << u_corr(q+(*kiter)) << " should be " << (vc(p) - vc0(p) + u(q+(*kiter))) << endl;
                                }
                            }
                            UNIT_TEST((abs(u_corr(q+(*kiter)) - (vc(p) - vc0(p) + u(q+(*kiter)))) < 1e-15));
                        }
                    }
                }

            }
        }
    }

    END_TEST();
    //==================================================================== 
    END_TEST_SUITE();
    #endif 
}
