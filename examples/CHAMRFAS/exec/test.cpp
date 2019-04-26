#include "LevelData.H"
#include "Multigrid.H"
#include "Proto.H"
#include "AMRFAS.H"

#define _LAPLACE_ 0
#define _MEHRSTELLEN_ 1

//#define OPERATOR _LAPLACE_
#define OPERATOR _MEHRSTELLEN_

#define GSRB TRUE

#include "BaseOp.H"
#include "LaplaceOp.H"
#include "MehrstellenOp.H"
#include "MehrstellenCorrectionOp.H"

using namespace std;

int main(int argc, char** argv)
{
#if CH_MPI
    MPI_Init(&argc, &argv);
#endif
    int TEST;
    int domainSize = 32;
    int numIter = 1;
    if (argc == 1)
    {
        cout << "Usage: " << endl; 
        cout << "\t Argument 1: Choose one of the following tests to run: " << endl;
        cout << "\t\tTest 0: Misc Testing (not guaranteed to do anything useful)" << endl;
        cout << "\t\tTest 1: AMR Multigrid" << endl;
        cout << "\t\tTest 2: AMR Operator" << endl;
        cout << "\t\tTest 3: Multigrid" << endl;
        cout << "\t\tTest 4: Boundary Interpolation" << endl;
        cout << "\t\tTest 5: Refluxing" << endl;
        cout << "\t Argument 2: Choose number of iterations (default 1). " << endl;
        cout << "\t Argument 3: Choose a coarse domainSize (default 32). " << endl;

        return 0;
    }
    if (argc >= 2) 
    {
        TEST = atoi(argv[1]);
    }
    if (argc >= 3) {
        numIter = atoi(argv[2]);
    }
    if (argc >= 4) {
        domainSize = atoi(argv[3]);
    }
    bool do_nesting = true;

#if OPERATOR==_LAPLACE_
    typedef LaplaceOp<FArrayBox> OP;
    std::cout << "Using Operator: LaplaceOp" << std::endl;
#elif OPERATOR==_MEHRSTELLEN_
    typedef MehrstellenOp<FArrayBox> OP;
    std::cout << "Using Operator: MehrstellenOp" << std::endl;
#endif
    typedef typename OP::patch BD;
    typedef FArrayBox DATA;

#ifdef GSRB
    std::cout << "GSRB = TRUE" << std::endl;
#else
    std::cout << "GSRB = FALSE" << std::endl;
#endif

     
    //====================================================================
    // Misc Testing
    if (TEST == 0)
    {
        cout << "Testing relaxation" << endl;
        Real L = 2.0*M_PI;
        Real dx = L/domainSize;
        auto domain = Proto::Box::Cube(domainSize);
        DisjointBoxLayout layout;
        buildLayout(layout, domain, Proto::Point::Ones());

        LevelData<DATA> Phi, R, Res;
        Phi.define(layout, OP::numcomps(), OP::ghost());
        R.define(layout, OP::numcomps(), Proto::Point::Zeros());
        Res.define(layout, OP::numcomps(), Proto::Point::Zeros());
        auto iter = layout.dataIterator();
        for (iter.begin(); iter.ok(); ++iter)
        {
            OP::patch rho = R[iter];
            Proto::forallInPlace_p(
                [=] PROTO_LAMBDA (Proto::Point& a_p, OP::var& a_data)
                {
                    Real x0 = a_p[0]*dx; Real x1 = x0 + dx;
                    a_data(0) = (sin(x1) - sin(x0))/dx;
                }, rho);
            OP::patch phi = Phi[iter];
            phi.setVal(0);
            OP::patch res = Res[iter];
            res.setVal(0);
        }

        OP op;
        op.define(layout, dx);
        writeLevel(Phi, "Test_relax_Phi_0.hdf5");
        writeLevel(R, "Test_relax_R_0.hdf5");
        for (int n = 0; n < numIter; n++) {
            op.relax(Phi, R, 1);
            
            op.residual(Res, Phi, R, dx);
            std::cout << "Residual: " << absMax(Res) << std::endl;
            writeLevel(Phi, "Test_relax_Phi_%i.hdf5", n+1);
            writeLevel(R, "Test_relax_R_%i.hdf5",n+1);
            writeLevel(Res, "Test_relax_Res_%i.hdf5",n+1);
        }
    
    }
    //====================================================================
    // AMRFAS Test
    if (TEST == 1)
    {
        cout << "Testing full AMRFAS algorithm" << endl;
        int numLevels = 2;
        Real L = 2.0*M_PI;
        Real cdx = L/domainSize;
        std::vector<DisjointBoxLayout> Layouts;
        
        std::vector<std::shared_ptr<LevelData<DATA>>> Phi;
        std::vector<std::shared_ptr<LevelData<DATA>>> Src;
        std::vector<std::shared_ptr<LevelData<DATA>>> Rhs;
        std::vector<std::shared_ptr<LevelData<DATA>>> Res;

        Layouts.resize(numLevels);
        Phi.resize(numLevels);
        Src.resize(numLevels);
        Rhs.resize(numLevels);
        Res.resize(numLevels);

#if OPERATOR==_MEHRSTELLEN_
        std::vector<std::shared_ptr<LevelData<DATA>>> RhsCorr;
        RhsCorr.resize(numLevels);
#endif

        Real dx = cdx;
        cout << "\tNumber of AMR Levels: " << numLevels << endl;
        cout << "\tReal size of domain: " << L << endl;
        cout << "\tNested domains: " << do_nesting << endl;
        cout << "\tTest Problem: Laplacian(<phi>) = <cos(x)>" << endl;
        cout << "\tInitial Phi: 0" << endl;
        for (int ii = 0; ii < numLevels; ii++)
        {
            auto& layout = Layouts[ii];
            if (do_nesting)
            {
                int s = ipow(AMR_REFRATIO,ii);
                Box domain = Proto::Box::Cube(domainSize/s);
                if (ii > 0)
                {
                    domain = domain.shift(Proto::Point::Ones(domainSize*0.5*(1-1.0/s)));
                }
                domain = domain.refine(s);
                std::cout << "\t\tDomain of level " << ii << ": " << domain << std::endl;
                if (ii == 0)
                {
                    buildLayout(layout, domain, Proto::Point::Ones());
                } else {
                    buildLayout(layout, domain, Proto::Point::Zeros());
                }
            } else {
                int s = ipow(AMR_REFRATIO,ii);
                Box domain = Proto::Box::Cube(domainSize);
                domain = domain.refine(s);
                std::cout << "\t\tDomain of level " << ii << ": " << domain << std::endl;
                buildLayout(layout, domain, Proto::Point::Ones());
            }
            std::cout << endl;
            Phi[ii] = std::make_shared<LevelData<DATA>>(layout, 1, OP::ghost());
            Src[ii] = std::make_shared<LevelData<DATA>>(layout, 1, OP::ghost());
#if OPERATOR==_LAPLACE_
            Rhs[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
#elif OPERATOR==_MEHRSTELLEN_
            Rhs[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Ones());
            RhsCorr[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
#endif            
            Res[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
            auto& phi = *(Phi[ii]);
            auto& src = *(Src[ii]);
            auto& rhs = *(Rhs[ii]);
            auto& res = *(Res[ii]);
            auto iter = phi.dataIterator();
            for (iter.reset(); iter.ok(); ++iter)
            {
                BD phi_i = phi[iter];
                BD rhs_i = rhs[iter];
                BD res_i = res[iter];
                res_i.setVal(0);
                forallInPlace_p(
                    [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_phi)
                    {
                        a_phi(0) = 0.0; 
                    }, phi_i);
                forallInPlace_p(
                    [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_rhs)
                    {
                        Real x0 = a_pt[0]*dx;
                        Real x1 = a_pt[0]*dx + dx;
                        a_rhs(0) = (sin(x1) - sin(x0))/dx; // = < cos(x) >
                    }, rhs_i);
            }
            phi.exchange();
            rhs.exchange();
            phi.copyTo(src);
            dx /= AMR_REFRATIO;
        } //end initialization
         
        AMRFAS<OP,DATA> amr_op(Layouts, dx*AMR_REFRATIO, numLevels-1, log2(1.0*domainSize) - 1);
        
        //do preprocessing
#if OPERATOR==_MEHRSTELLEN_
        AMRFAS<MehrstellenCorrectionOp<DATA>,DATA> correctOp(Layouts, dx*AMR_REFRATIO, numLevels-1, 1);
        correctOp(RhsCorr, Rhs);
        std::cout << "Integral of correction: " << integrate(RhsCorr, cdx) << endl;
        std::cout << "Integral of non-conditioned RHS: " << integrate(Rhs, cdx) << endl;
        for (int ii = 0; ii < numLevels; ii++)
        {
            LevelData<DATA>& rhs = *(Rhs[ii]);
            LevelData<DATA>& corr = *(RhsCorr[ii]);
            auto iter = rhs.dataIterator();
            for (iter.begin(); iter.ok(); ++iter)
            {
                OP::patch rhsPatch = rhs[iter];
                OP::patch corrPatch = corr[iter];
                corrPatch += rhsPatch;
            }
        }
        correctOp.write(RhsCorr, "AMR_RhsCorrected.hdf5");
        std::cout << "Integral of conditioned RHS: " << integrate(RhsCorr, cdx) << endl;
#else
        cout << "Integral of RHS: " << integrate(Rhs, cdx) << endl; 
#endif       
        amr_op.write(Src, "AMR_Src.hdf5");
        amr_op.write(Rhs, "AMR_Rhs.hdf5");
        cout << "Integral of Res: " << integrate(Res, cdx) << endl; 
        cout << "Integral of Phi: " << integrate(Phi, cdx) << endl; 
        for (int nn = 0; nn < numIter; nn++)
        {
            amr_op.write(Res, "AMR_Res.%i.hdf5", nn);
            amr_op.write(Phi, "AMR_Phi.%i.hdf5", nn);
#if OPERATOR==_MEHRSTELLEN_
            amr_op.vcycle(Phi, RhsCorr, Res);
#else 
            amr_op.vcycle(Phi, Rhs, Res);
#endif
            cout << scientific << "Residual: Max = " << absMax(Res);
            cout << "\t\tIntegral(Res) = " << integrate(Res, cdx) << endl;
        }
        amr_op.write(Res, "AMR_Res.%i.hdf5", numIter);
        amr_op.write(Phi, "AMR_Phi.%i.hdf5", numIter);

        for (int level = 0; level < numLevels; ++level)
        {
            auto& rhs = *(Rhs[level]);
            auto& src = *(Src[level]);
            auto& phi = *(Phi[level]);
            auto iter = src.dataIterator();
            for (iter.reset(); iter.ok(); ++iter)
            {
                BD srcPatch = src[iter];
                BD phiPatch = phi[iter];
                BD rhsPatch = rhs[iter];
                rhsPatch += phiPatch; //because the solution for phi is -<cos(x)> = -rhs
                phiPatch -= srcPatch;
            }
        }
        amr_op.write(Rhs, "AMR_Error.hdf5");
        amr_op.write(Phi, "AMR_Diff.hdf5");
        Real innerError = absMax(Rhs, 0, true, 1);
        Real bdryError = absMax(Rhs, 0, false, 1);

        cout << "Interior Error: " << innerError << endl;
        cout << "Boundary Error: " << bdryError << endl << endl;

    } // End AMRFAS test 
    //====================================================================
    // AMR Operator Test
    if (TEST == 2)
    {
        cout << "Testing Simple AMROperator" << endl;
        const int numLevels = 2;
        Real ei[numIter];
        Real eb[numIter];
        cout << "\tNumber of AMR Levels: " << numLevels << endl;
        cout << "\tNested Domains: " << true << endl; // no other option at the moment
        for (int n = 0; n < numIter; n++)
        {
            Real cdx = 2.0*M_PI/domainSize;
            std::vector<DisjointBoxLayout> Layouts;
            std::vector<std::shared_ptr<LevelData<DATA>>> Phi;
            std::vector<std::shared_ptr<LevelData<DATA>>> LPhi;
            std::vector<std::shared_ptr<LevelData<DATA>>> Rhs;
#if OPERATOR==_MEHRSTELLEN_
            std::vector<std::shared_ptr<LevelData<DATA>>> RhsCorr;
#endif
            Layouts.resize(numLevels);
            Phi.resize(numLevels);
            LPhi.resize(numLevels);
            Rhs.resize(numLevels);
#if OPERATOR==_MEHRSTELLEN_
            RhsCorr.resize(numLevels);
#endif
            Real dx = cdx;
            for (int ii = 0; ii < numLevels; ii++)
            {
                auto& layout = Layouts[ii];
                if (do_nesting)
                {
                    int s = ipow(AMR_REFRATIO,ii);
                    Box domain = Proto::Box::Cube(domainSize/s);
                    if (ii > 0)
                    {
                        domain = domain.shift(Proto::Point::Ones(domainSize*0.5*(1-1.0/s)));
                    }
                    domain = domain.refine(s);
                    std::cout << "\t\tDomain of level " << ii << ": " << domain << std::endl;
                    if (ii == 0)
                    {
                        buildLayout(layout, domain, Proto::Point::Ones());
                    } else {
                        buildLayout(layout, domain, Proto::Point::Zeros());
                    }
                } else {
                    int s = ipow(AMR_REFRATIO,ii);
                    Box domain = Proto::Box::Cube(domainSize);
                    domain = domain.refine(s);
                    std::cout << "\t\tDomain of level " << ii << ": " << domain << std::endl;
                    buildLayout(layout, domain, Proto::Point::Ones());
                }
                Phi[ii] = std::make_shared<LevelData<DATA>>(layout, 1, OP::ghost());
                LPhi[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
                Rhs[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Ones());
#if OPERATOR==_MEHRSTELLEN_
                RhsCorr[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
#endif
                auto& phi = *(Phi[ii]);
                auto& Lphi = *(LPhi[ii]);
                auto& rhs = *(Rhs[ii]);
                auto iter = phi.dataIterator();
                for (iter.reset(); iter.ok(); ++iter)
                {
                    BD phi_i = phi[iter];
                    BD Lphi_i = Lphi[iter];
                    BD rhs_i = rhs[iter];
                    forallInPlace_p(
                            [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_phi, OP::var& a_rhs)
                            {
                                Real x0 = a_pt[0]*dx;
                                Real x1 = x0 + dx;
                                Real y0 = a_pt[1]*dx;
                                Real y1 = y0 + dx;
                                a_phi(0) = (sin(x1) - sin(x0))/dx;//*(sin(y1) - sin(y0))/(dx*dx);
                                a_rhs(0) = -a_phi(0);
                            }, phi_i, rhs_i);
                    Lphi_i.setVal(0);
                }
                phi.exchange();
                rhs.exchange();
                dx /= AMR_REFRATIO;
            } //end initialization
            //do preprocessing
#if OPERATOR==_MEHRSTELLEN_
            // build an AMR Laplace operator to help to the preprocessing
            AMRFAS<MehrstellenCorrectionOp<DATA>,DATA> correctOp(Layouts, dx*AMR_REFRATIO, numLevels-1, 1);
            correctOp(RhsCorr, Rhs);
            std::cout << "Integral of correction: " << integrate(RhsCorr, cdx) << endl;
            std::cout << "Integral of non-conditioned RHS: " << integrate(Rhs, cdx) << endl;
            for (int ii = 0; ii < numLevels; ii++)
            {
                LevelData<DATA>& rhs = *(Rhs[ii]);
                LevelData<DATA>& corr = *(RhsCorr[ii]);
                auto iter = rhs.dataIterator();
                for (iter.begin(); iter.ok(); ++iter)
                {
                    OP::patch rhsPatch = rhs[iter];
                    OP::patch corrPatch = corr[iter];
                    corrPatch += rhsPatch;
                }
            }
            correctOp.write(RhsCorr, "AMR_RhsCorrected.%i.hdf5",n);
            std::cout << "Integral of conditioned RHS: " << integrate(RhsCorr, cdx) << endl;
#endif       

            AMRFAS<OP,DATA> amr_op(Layouts, dx*AMR_REFRATIO , numLevels-1, 1);

            cout << "Integral of Phi: " << integrate(Phi, cdx) << endl;
            amr_op(LPhi, Phi);
            
            cout << "Integral of L(Phi): " << integrate(LPhi, cdx) << endl; 
            amr_op.write(LPhi, "AMR_LPhi.%i.hdf5",n);
            amr_op.write(Phi, "AMR_Phi.%i.hdf5",n);
            for (int level = 0; level < numLevels; ++level)
            {
#if OPERATOR==_MEHRSTELLEN_
                auto& rhs = *(RhsCorr[level]);
                std::cout << "Computing error using corrected Rho" << std::endl;
#else
                auto& rhs = *(Rhs[level]);
#endif
                auto& Lphi = *(LPhi[level]);
                auto iter = Lphi.dataIterator();
                for (iter.reset(); iter.ok(); ++iter)
                {
                    BD LphiPatch = Lphi[iter];
                    BD rhsPatch = rhs[iter];
                    LphiPatch -= rhsPatch;
                }
            }
            amr_op.write(LPhi, "AMR_Error.%i.hdf5",n);
            Real innerError = absMax(LPhi, 0, true, 1);
            Real bdryError = absMax(LPhi, 0, false, 1);
            ei[n] = innerError;
            eb[n] = bdryError;

            cout << "Interior Error: " << innerError << endl;
            cout << "Boundary Error: " << bdryError << endl << endl;
            domainSize *= 2;
        }
        
        for (int ii = 1; ii < numIter; ii++)
        {
            Real ri = log2(ei[ii-1]/ei[ii]);
            Real rb = log2(eb[ii-1]/eb[ii]);
            cout << "Interior Rate: " << ri << endl;
            cout << "Boundary Rate: " << rb << endl << endl;
        }
    } // End AMRFAS test 
    //====================================================================
    //====================================================================
    else if (TEST == 3)
    {
        bool doAMR = false;
        int numLevels = 2;
        //int numLevels = log(domainSize*1.0)/log(2.0)-1;
        Real cdx = 2.0*M_PI/domainSize;
        Real dx = cdx / 2.0; 

        std::cout << "Running Multigrid:" << std::endl;
        std::cout << "\tTesting AMR: " << doAMR << std::endl;
        std::cout << "\tDomain Size: " << domainSize << std::endl;
        std::cout << "\tMax Box Size: " << MAXBOXSIZE << std::endl;
        std::cout << "\tNumber of Multigrid Levels: " << numLevels << std::endl;
            
        auto coarseDomain = Proto::Box::Cube(domainSize);
        auto coarseFineDomain = Proto::Box::Cube(domainSize/2).shift(Proto::Point::Ones(domainSize/4));
        auto fineDomain = coarseFineDomain.refine(AMR_REFRATIO);
        
        cout << "Coarse Domain: " << coarseDomain << endl;
        cout<< "Coarsened Fine Domain: " << coarseFineDomain << endl;
        cout<< "Fine Domain: " << fineDomain << endl;

        DisjointBoxLayout fineLayout, coarseLayout, coarseTemp;
        buildLayout(fineLayout, fineDomain, Proto::Point::Zeros());
        buildLayout(coarseLayout, coarseDomain, Proto::Point::Ones());
        coarsen_dbl(coarseTemp, fineLayout, AMR_REFRATIO);

        LevelData<DATA> PhiC(coarseLayout, OP::numcomps(), IntVect::Unit);
        LevelData<DATA> RC(coarseLayout, OP::numcomps(), IntVect::Zero);
        LevelData<DATA> ResC(coarseLayout, OP::numcomps(), IntVect::Zero);
        LevelData<DATA> Phi(fineLayout, OP::numcomps(), IntVect::Unit);
        LevelData<DATA> R(fineLayout, OP::numcomps(), IntVect::Zero);
        LevelData<DATA> Res(fineLayout, OP::numcomps(), IntVect::Zero);
    
        auto fiter = fineLayout.dataIterator();
        for (fiter.begin(); fiter.ok(); ++fiter)
        {
            BD res = Res[fiter];
            res.setVal(0);
            BD phi = Phi[fiter];
            phi.setVal(0);
            BD rhs = R[fiter];
            Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, Proto::Var<Real>& a_rhs)
            {
                Real x0 = a_p[0]*dx;
                Real x1 = x0 + dx;
                a_rhs(0) = (sin(x1) - sin(x0))/dx;
            }, rhs);
        }

        auto citer = coarseLayout.dataIterator();
        for (citer.begin(); citer.ok(); ++citer)
        {
            BD resC = ResC[citer];
            resC.setVal(0);
            BD phiC = PhiC[citer];
            phiC.setVal(0);
            BD rhsC = RC[citer];
            Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, Proto::Var<Real>& a_rhs)
            {
                Real x0 = a_p[0]*cdx;
                Real x1 = x0 + cdx;
                a_rhs(0) = (sin(x1) - sin(x0))/cdx;
            }, rhsC);
        }

        if (doAMR)
        {
            std::cout << "AMR test is not fully functional" << std::endl;
            Multigrid<OP, DATA> amr_mg(fineLayout, dx, AMR_REFRATIO/MG_REFRATIO - 1, true, 1);
            amr_mg.vcycle(Phi,PhiC,R);
        } else {
            Multigrid<OP, DATA> mg(coarseLayout, cdx, numLevels-1, false);

            double resnorm = 0.0;
            OP op;
            op.define(coarseLayout,cdx);
            for (int ii = 0; ii < numIter; ii++)
            {
                writeLevel(PhiC, "MG_Phi.%i.hdf5",ii);
                mg.vcycle(PhiC,RC); 
                resnorm = op.residual(ResC,PhiC,RC);
                std::cout << scientific << "iteration number = " << ii << ", Residual norm: " << resnorm << std::endl;
            }
        }
    } // End Multigrid test
    else if (TEST == 4) {
        Real L = 2.0*M_PI;

        Real error[numIter]; 
        for (int nn = 0; nn < numIter; nn++)
        {
            Real cdx = L/domainSize;
            Real dx = cdx/AMR_REFRATIO;
            std::cout << "Testing Boundary Interpolation" << std::endl;
            PROTO_ASSERT(domainSize % 4 == 0, "domainSize should be a multiple of 4");
            auto coarseDomain = Proto::Box::Cube(domainSize);
            auto coarseFineDomain = Proto::Box::Cube(domainSize/2).shift(Proto::Point::Ones(domainSize/4));
            auto fineDomain = coarseFineDomain.refine(AMR_REFRATIO);

            DisjointBoxLayout fineLayout, coarseFineLayout, coarseLayout;
            buildLayout(coarseLayout, coarseDomain, Proto::Point::Ones());
            buildLayout(fineLayout, fineDomain, Proto::Point::Zeros());
            coarsen_dbl(coarseFineLayout, fineLayout, AMR_REFRATIO);   

            LevelData<DATA> PhiC(coarseLayout, OP::numcomps(), OP::ghost());
            LevelData<DATA> Phi(fineLayout, OP::numcomps(), OP::ghost());
            LevelData<DATA> Soln(fineLayout, OP::numcomps(), OP::ghost());        
            auto citer = coarseLayout.dataIterator();
            for (citer.begin(); citer.ok(); ++citer)
            {
                OP::patch phiC = PhiC[citer];
                Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, OP::var& a_data)
                        {
                        Real x0 = a_p[0]*cdx;
                        Real x1 = x0 + cdx;
                        a_data(0) = (sin(x1) - sin(x0))/cdx;
                        }, phiC);
            }
            auto fiter = fineLayout.dataIterator();
            for (fiter.begin(); fiter.ok(); ++fiter)
            {
                OP::patch phi = Phi[fiter];
                phi.setVal(0);
                OP::patch soln = Soln[fiter];
                Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, OP::var& a_data)
                        {
                        if (!fineDomain.contains(a_p))
                        {
                        Real x0 = a_p[0]*dx;
                        Real x1 = x0 + dx;
                        a_data(0) = (sin(x1) - sin(x0))/dx;
                        } else {
                        a_data(0) = 0;
                        }
                        }, soln);
            }

            OP op;
            op.define(fineLayout, dx, true); 
            op.interpBoundary(Phi, PhiC);       

            writeLevel(PhiC, "InterpTest_PhiC.%i.hdf5",nn);
            writeLevel(Phi, "InterpTest_Phi.%i.hdf5",nn);
            writeLevel(Soln, "InterpTest_Soln.%i.hdf5",nn);

            error[nn] = 0.0;
            for (fiter.reset(); fiter.ok(); ++fiter)
            {
                OP::patch soln = Soln[fiter];
                OP::patch phi = Phi[fiter];
                soln -= phi;
                error[nn] = std::max(error[nn], soln.absMax());
            }
            std::cout << "Error: " << error[nn] << std::endl;
            domainSize *= 2;
        }
        for (int ii = 1; ii < numIter; ii++)
        {
            std::cout << "Rate: " << log2(error[ii-1]/error[ii]) << std::endl;
        }
    } // End of Interpolation test
    else if (TEST == 5) {
        std::cout << "Testing Refluxing" << std::endl;
        Real L = 2.0*M_PI;

        Real cdx = L/domainSize;
        Real dx = cdx/AMR_REFRATIO;
        
        PROTO_ASSERT(domainSize % 4 == 0, "domainSize should be a multiple of 4");
        auto coarseDomain = Proto::Box::Cube(domainSize);
        auto coarseFineDomain = Proto::Box::Cube(domainSize/2).shift(Proto::Point::Ones(domainSize/4));
        auto fineDomain = coarseFineDomain.refine(AMR_REFRATIO);

        DisjointBoxLayout fineLayout, coarseFineLayout, coarseLayout;
        buildLayout(coarseLayout, coarseDomain, Proto::Point::Ones());
        buildLayout(fineLayout, fineDomain, Proto::Point::Zeros());
        coarsen_dbl(coarseFineLayout, fineLayout, AMR_REFRATIO);   

        LevelData<DATA> PhiC(coarseLayout, OP::numcomps(), OP::ghost());
        LevelData<DATA> LPhiC(coarseLayout, OP::numcomps(), Proto::Point::Zeros());
        LevelData<DATA> NullPhi(coarseLayout, OP::numcomps(), OP::ghost());
        LevelData<DATA> Phi(fineLayout, OP::numcomps(), OP::ghost());
        LevelData<DATA> LPhi(fineLayout, OP::numcomps(), Proto::Point::Zeros());
        LevelData<DATA> Flux(coarseLayout, OP::numcomps(), Proto::Point::Zeros());        
        auto citer = coarseLayout.dataIterator();
        for (citer.begin(); citer.ok(); ++citer)
        {
            OP::patch phiC = PhiC[citer];
            Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, OP::var& a_data)
                    {
                    Real x0 = a_p[0]*cdx;
                    Real x1 = x0 + cdx;
                    a_data(0) = (sin(x1) - sin(x0))/cdx;
                    }, phiC);
            OP::patch flux = Flux[citer];
            OP::patch nullPhi = NullPhi[citer];
            flux.setVal(0);
            nullPhi.setVal(0);
        }
        auto fiter = fineLayout.dataIterator();
        for (fiter.begin(); fiter.ok(); ++fiter)
        {
            OP::patch phi = Phi[fiter];
            Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, OP::var& a_data)
                    {
                    Real x0 = a_p[0]*dx;
                    Real x1 = x0 + dx;
                    a_data(0) = (cos(x1) - cos(x0))/dx;
                    }, phi);
        }
        
        LevelFluxRegister Register;
        Register.define(fineLayout, coarseLayout, fineLayout.physDomain(), AMR_REFRATIO, OP::numcomps());

        OP op;
        op.define(fineLayout, dx, true); 
        op.reflux(Flux, PhiC, Phi, Register);
        op.apply(LPhi, Phi, dx);
        op.apply(LPhiC, PhiC, cdx);
        // reset interior garbage
        for (citer.reset(); citer.ok(); ++citer)
        {
            OP::patch flux = Flux[citer];
            OP::patch LphiC = LPhiC[citer];
            Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, OP::var& a_flux, OP::var& a_phiC)
            {
                if (coarseFineDomain.contains(a_p))
                {
                    a_flux(0) = 0;
                    a_phiC(0) = 0;
                }
            }, flux, LphiC);
        }

        writeLevel(PhiC, "RefluxTest_PhiC.hdf5");
        writeLevel(Phi, "RefluxTest_Phi.hdf5");
        writeLevel(Flux, "RefluxTest_Soln.hdf5");

        Real intFlux = integrate(Flux, cdx);
        Real intCoarse = integrate(LPhiC, cdx);
        Real intFine = integrate(LPhi, dx);

        std::cout << "Integral of L(Phi) on fine grid: " << intFine << std::endl;
        std::cout << "Integral of L(PhiC) on coarse grid: " << intCoarse << std::endl;
        std::cout << "Sum: " << intFine + intCoarse << std::endl;
        std::cout << "Integral of composite register: " << intFlux << std::endl;
        std::cout << "Error: " << std::abs(intFine + intCoarse - intFlux) << std::endl;
    }
#if CH_MPI
    CH_TIMER_REPORT();
    MPI_Finalize();
#endif
}
