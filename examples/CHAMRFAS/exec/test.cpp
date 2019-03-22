#include "LevelData.H"
#include "Multigrid.H"
#include "Proto.H"
#include "AMRFAS.H"

#include "BaseOp.H"
#include "TestOp.H"

//using namespace Proto;
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
        cout << "\t\tTest 4: Refactored Operator" << endl;
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

    typedef TestOp<FArrayBox> OP;
    typedef typename OP::patch BD;
    typedef FArrayBox DATA;

     
    //====================================================================
    // Misc Testing
    if (TEST == 0)
    {
        cout << "Currently Testing: GSRB" << endl;
        
        Real L = 2.0*M_PI;
        int numIter = 1;
        Real error[numIter];
        for (int n = 0; n < numIter; n++)
        {
            cout << "Iteration: " << n << endl;
            Real cdx = L/domainSize;
            Real dx = cdx/AMR_REFRATIO;

            auto coarseDomain = Proto::Box::Cube(domainSize);
            auto coarseFineDomain = Proto::Box::Cube(domainSize/2).shift(Proto::Point::Ones(domainSize/4));
            auto fineDomain = coarseFineDomain.refine(AMR_REFRATIO);
            cout << "Coarse Domain: " << coarseDomain << endl;
            cout<< "Coarsened Fine Domain: " << coarseFineDomain << endl;
            cout<< "Fine Domain: " << fineDomain << endl;

            DisjointBoxLayout fineLayout, coarseLayout, coarseTemp;
            buildLayout(fineLayout, fineDomain, Proto::Point::Zeros());
            buildLayout(coarseLayout, coarseDomain, Proto::Point::Ones());
            coarsen(coarseTemp, fineLayout, AMR_REFRATIO);
            auto fiter = fineLayout.dataIterator();
            auto citer = coarseLayout.dataIterator();
        
            LevelData<DATA> Phi(coarseLayout, 1, Proto::Point::Ones());
            LevelData<DATA> Rhs(coarseLayout, 1, Proto::Point::Zeros());
            LevelData<DATA> Res(coarseLayout, 1, Proto::Point::Zeros());
            for (citer.begin(); citer.ok(); ++citer)
            {
                BD phi = Phi[citer];
                BD rhs = Rhs[citer];
                BD res = Res[citer];
                res.setVal(0);
                forallInPlace_p(
                    [=] PROTO_LAMBDA (Proto::Point& a_pt, Proto::Var<Real>& a_src, Proto::Var<Real>& a_rhs)
                    {
                        Real x0 = a_pt[0]*cdx;
                        Real x1 = a_pt[0]*cdx + cdx;
                        a_rhs(0) = (sin(x1) - sin(x0))/cdx; // = < cos(x) >
                        a_src(0) = (-1.0)*a_rhs(0) + 0.1*(sin(5*x1) - sin(5*x0))/cdx;
                    }, phi, rhs
                );
            }
            
            Phi.exchange();
            
            TestOp<DATA> op;
            op.define(coarseLayout, cdx, false);
            int N = 10;
            for (int ii = 0; ii < N; ii++)
            {
                writeLevel(Phi, "GSRB_TEST_Phi.%i.hdf5", ii);
                op.residual(Res, Phi, Rhs);
                cout << "Max of Residual: " << absMax(Res) << endl;
                op.relax(Phi, Rhs, 1);
            }
            writeLevel(Phi, "GSRB_TEST_Phi.%i.hdf5",N);
             
            domainSize *= 2;
        }
        for (int ii = 1; ii < numIter; ii++)
        {
            cout << "Convergence Rate: " << log2(error[ii-1]/error[ii]) << endl;
        }
    }

    //====================================================================
    // AMRFAS Test
    if (TEST == 1)
    {
        cout << "Testing full AMRFAS algorithm" << endl;
        int numLevels = 2;
        bool do_nesting = true;
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
            Phi[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Ones(1));
            Src[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Ones(1));
            Rhs[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
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
            phi.copyTo(src);
            dx /= AMR_REFRATIO;
        } //end initialization
        
        AMRFAS<OP,DATA> amr_op(Layouts, dx*AMR_REFRATIO , numLevels-1, log2(1.0*domainSize) - 1);
        amr_op.write(Src, "AMR_Src.hdf5");
        amr_op.write(Rhs, "AMR_Rhs.hdf5");
        cout << "Integral of RHS: " << integrate(Rhs, cdx) << endl; 
        cout << "Integral of Res: " << integrate(Res, cdx) << endl; 
        cout << "Integral of Phi: " << integrate(Phi, cdx) << endl; 
        for (int nn = 0; nn < numIter; nn++)
        {
            amr_op.write(Res, "AMR_Res.%i.hdf5", nn);
            amr_op.write(Phi, "AMR_Phi.%i.hdf5", nn);
            amr_op.vcycle(Phi, Rhs, Res);
            cout << scientific << "Residual: Max = " << absMax(Res) << "\t\tIntegral = " << integrate(Res, cdx) << endl;
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
        Real innerError = absMax(Rhs, true, 1);
        Real bdryError = absMax(Rhs, false, 1);

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
            std::vector<std::shared_ptr<LevelData<DATA>>> Src;
            std::vector<std::shared_ptr<LevelData<DATA>>> Rhs;

            Layouts.resize(numLevels);
            Phi.resize(numLevels);
            Src.resize(numLevels);
            Rhs.resize(numLevels);
            Real dx = cdx;
            for (int ii = 0; ii < numLevels; ii++)
            {
                int s = ipow(AMR_REFRATIO,ii);
                Box domain = Proto::Box::Cube(domainSize/s);

                if (ii > 0)
                {
                    domain = domain.shift(Proto::Point::Ones(domainSize*0.5*(1-1.0/s)));
                }

                domain = domain.refine(s);
                auto& layout = Layouts[ii];
                if (ii == 0)
                {
                    buildLayout(layout, domain, Proto::Point::Ones());
                } else {
                    buildLayout(layout, domain, Proto::Point::Zeros());
                }
                Phi[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Ones(1));
                Src[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Ones(1));
                Rhs[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
                auto& phi = *(Phi[ii]);
                auto& src = *(Src[ii]);
                auto& rhs = *(Rhs[ii]);
                auto iter = phi.dataIterator();
                for (iter.reset(); iter.ok(); ++iter)
                {
                    BD phi_i = phi[iter];
                    BD rhs_i = rhs[iter];
                    forallInPlace_p(
                            [=] PROTO_LAMBDA (Proto::Point& a_pt, Proto::Var<Real>& a_phi)
                            {
                                a_phi(0) = (-cos(a_pt[0]*dx+dx) + cos(a_pt[0]*dx))/dx;
                            }, phi_i);
                    rhs_i.setVal(0);
                }
                phi.exchange();
                phi.copyTo(src);
                dx /= AMR_REFRATIO;
            } //end initialization

            AMRFAS<OP,DATA> amr_op(Layouts, dx*AMR_REFRATIO , numLevels-1, 1);
            amr_op(Rhs, Phi);
            cout << "Integral of L(Phi): " << integrate(Rhs, cdx) << endl; 
            if (n == 0)
            {
                amr_op.write(Rhs, "AMR_Rhs.hdf5");
                amr_op.write(Phi, "AMR_Phi.hdf5");
                amr_op.write(Src, "AMR_Src.hdf5");
            }
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
                    rhsPatch += srcPatch;
                    phiPatch -= srcPatch;
                }
            }
            if (n == 0)
            {
                amr_op.write(Rhs, "AMR_Error.hdf5");
                amr_op.write(Phi, "AMR_Diff.hdf5");
            }
            Real innerError = absMax(Rhs, true, 1);
            Real bdryError = absMax(Rhs, false, 1);
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
        int numLevels = log(domainSize*1.0)/log(2.0)-1;
        Real cdx = 2.0*M_PI/domainSize;
        Real dx = cdx / 2.0; 

        std::cout << "Running Multigrid:" << std::endl;
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
        Real dx = L/domainSize;
        std::cout << "Testing refactored operator code" << std::endl;
        auto domain = Proto::Box::Cube(domainSize);
        DisjointBoxLayout layout;
        buildLayout(layout, domain, Proto::Point::Ones());
        
        LevelData<DATA> Phi(layout, 1, Proto::Point::Ones());
        LevelData<DATA> Rhs(layout, 1, Proto::Point::Zeros());
        
        auto iter = layout.dataIterator();
        for (iter.begin(); iter.ok(); ++iter)
        {
            BD phi = Phi[iter];
            phi.setVal(0);
            BD rhs = Rhs[iter];
            Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, Proto::Var<Real>& a_data)
            {
                Real x0 = a_p[0]*dx;
                Real x1 = x0 + dx;
                a_data(0) = (sin(x1) - sin(x0))/dx;
            }, rhs);
        }

        Multigrid<OP, DATA> mg(layout, dx, 1, false);

    }
#if CH_MPI
    CH_TIMER_REPORT();
    MPI_Finalize();
#endif
}
