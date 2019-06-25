#include "LevelData.H"
#include "Multigrid.H"
#include "Proto.H"
#include "AMRFAS.H"

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
        cout << "\t\tTest 6: Residual" << endl;
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

    std::cout << "Coarse Domain Size: " << domainSize << std::endl;
    std::cout << "AMR ref ratio: " << AMR_REFRATIO << std::endl;
    std::cout << "MG ref ratio: " << MG_REFRATIO << std::endl;
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
#if OPERATOR==_MEHRSTELLEN_
        LevelData<DATA> RCorr;
        R.define(layout, OP::numcomps(), Proto::Point::Ones());
        RCorr.define(layout, OP::numcomps(), Proto::Point::Zeros());
#endif
        Res.define(layout, OP::numcomps(), Proto::Point::Zeros());
        auto iter = layout.dataIterator();
        for (iter.begin(); iter.ok(); ++iter)
        {
            OP::patch rho = R[iter];
            OP::patch phi = Phi[iter];
            OP::patch res = Res[iter];
            res.setVal(0);
            Proto::forallInPlace_p(
                [=] PROTO_LAMBDA (Proto::Point& a_p, OP::var& a_data)
                {
                    Real x0 = a_p[0]*dx; Real x1 = x0 + dx;
                    a_data(0) = (sin(x1) - sin(x0))/dx;
                }, rho);
            Proto::forallInPlace_p(
                [=] PROTO_LAMBDA (Proto::Point& a_p, OP::var& a_data)
                {
                    Real x0 = a_p[0]*dx; Real x1 = x0 + dx;
                    a_data(0) = -(sin(x1) - sin(x0))/dx;
                }, phi);
#if OPERATOR==_MEHRSTELLEN_
            auto S = 1.0*Proto::Shift::Zeros() + 1.0/12.0*Proto::Stencil<Real>::Laplacian();
            OP::patch corr = RCorr[iter];
            corr |= S(rho);
#endif
        }

        OP op;
        op.define(layout, dx);
#if OPERATOR==_MEHRSTELLEN_
        op.residual(Res, Phi, RCorr, dx);
#else
        op.residual(Res, Phi, R, dx);
#endif
        writeLevel(Res, "TestResidual.hdf5");
        std::cout << "Residual: " << absMax(Res) << std::endl;
    
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
        std::vector<std::shared_ptr<LevelData<DATA>>> Sln;
        std::vector<std::shared_ptr<LevelData<DATA>>> Rhs;
        std::vector<std::shared_ptr<LevelData<DATA>>> Res;
        std::vector<std::shared_ptr<LevelData<DATA>>> Err;

        Layouts.resize(numLevels);
        Phi.resize(numLevels);
        Sln.resize(numLevels);
        Rhs.resize(numLevels);
        Res.resize(numLevels);
        Err.resize(numLevels);

#if OPERATOR==_MEHRSTELLEN_
        std::vector<std::shared_ptr<LevelData<DATA>>> RhsCorr;
        RhsCorr.resize(numLevels);
#endif

        Real dx = cdx;
        cout << "\tNumber of AMR Levels: " << numLevels << endl;
        cout << "\tReal size of domain: " << L << endl;
        cout << "\tTest Problem: Laplacian(<phi>) = <cos(x)>" << endl;
        cout << "\tInitial Phi: 0" << endl;
        for (int ii = 0; ii < numLevels; ii++)
        {
            auto& layout = Layouts[ii];
            int s = ipow(AMR_REFRATIO,ii);
#if NESTING==_FULL_NESTING_
            Box domain = Proto::Box::Cube(domainSize/s);
            std::cout << "Nesting in all directions" << std::endl;
            std::cout << "\t\tUnshifted Domain: " << domain << std::endl;
            if (ii > 0)
            {
                domain = domain.shift(Proto::Point::Ones(domainSize*0.5*(1-1.0/s)));
            }
            std::cout << "\t\tShifted Domain: " << domain << std::endl;
            domain = domain.refine(s);
            std::cout << "\t\tShifted and Refined Domain: " << domain << std::endl;
            std::cout << "Building layout..." << std::endl;
            if (ii == 0)
            {
                buildLayout(layout, domain, Proto::Point::Ones());
            } else {
                buildLayout(layout, domain, Proto::Point::Zeros());
            }
#elif NESTING==_X_NESTING_
            Box domain = Proto::Box(Proto::Point(domainSize/s-1, domainSize-1));
            std::cout << "Nesting X Direction" << std::endl;
            std::cout << "\t\tUnshifted Domain: " << domain << std::endl;
            if (ii > 0)
            {
                domain = domain.shift(Proto::Point::Basis(0, domainSize*0.5*(1-1.0/s)));
            }
            std::cout << "\t\tShifted Domain: " << domain << std::endl;
            domain = domain.refine(s);
            std::cout << "\t\tShifted and Refined Domain: " << domain << std::endl;
            std::cout << "Building layout..." << std::endl;
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
            std::cout << "\t\tDomain of level " << ii << ": " << domain << std::endl;
            std::cout << endl;
            Phi[ii] = std::make_shared<LevelData<DATA>>(layout, 1, OP::ghost());
            Sln[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
            Err[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
            Res[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
#if OPERATOR==_LAPLACE_
            Rhs[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
#elif OPERATOR==_MEHRSTELLEN_
            Rhs[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Ones());
            RhsCorr[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
#endif            
            auto& phi = *(Phi[ii]);
            auto& sln = *(Sln[ii]);
            auto& rhs = *(Rhs[ii]);
            auto& res = *(Res[ii]);
            auto& err = *(Err[ii]);
            auto iter = phi.dataIterator();
            for (iter.reset(); iter.ok(); ++iter)
            {
                BD phi_i = phi[iter];
                BD rhs_i = rhs[iter];
                BD res_i = res[iter];
                BD sln_i = sln[iter];
                BD err_i = err[iter];
                phi_i.setVal(0);
                res_i.setVal(0);
                err_i.setVal(0);
                forallInPlace_p(
                    [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_rhs)
                    {
                        Real x0 = a_pt[0]*dx;
                        Real x1 = a_pt[0]*dx + dx;
                        a_rhs(0) = (sin(x1) - sin(x0))/dx; // = < cos(x) >
                    }, rhs_i);
                forallInPlace_p(
                    [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_sln)
                    {
                        Real x0 = a_pt[0]*dx;
                        Real x1 = a_pt[0]*dx + dx;
                        a_sln(0) = -(sin(x1) - sin(x0))/dx; // = < -cos(x) >
                    }, sln_i);
                /*
                forallInPlace_p(
                    [=] PROTO_LAMBDA (Proto::Point& a_pt, OP::var& a_phi)
                    {
                        Real x0 = a_pt[0]*dx;
                        Real x1 = a_pt[0]*dx + dx;
                        a_phi(0) = -(sin(x1) - sin(x0))/dx; // = < -cos(x) >
                    }, phi_i);
                */
            }
            phi.exchange();
            rhs.exchange();
            dx /= AMR_REFRATIO;
        } //end initialization
         
        AMRFAS<OP,DATA> amr_op(Layouts, dx*AMR_REFRATIO, log2(1.0*domainSize) - 1);
        
#if OPERATOR==_MEHRSTELLEN_
        AMRFAS<MehrstellenCorrectionOp<DATA>,DATA> correctOp(Layouts, dx*AMR_REFRATIO, 1);
        correctOp(RhsCorr, Rhs);
        correctOp.write(RhsCorr, "AMR_RhoCorrection.hdf5");
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
        cout << "Integral of Phi: " << integrate(Phi, cdx) << endl; 
        amr_op.write(Sln, "AMR_Sln.hdf5");
        amr_op.write(Rhs, "AMR_Rhs.hdf5");
#if OPERATOR==_MEHRSTELLEN_
        amr_op.residual(Res, Phi, RhsCorr); 
#else
        amr_op.residual(Res, Phi, Rhs); 
#endif
        cout << "Initial Integral(Res) = " << integrate(Res, cdx) << endl;
        for (int nn = 0; nn < numIter; nn++)
        {
            amr_op.write(Res, "AMR_Res.%i.hdf5", nn);
            amr_op.write(Phi, "AMR_Phi.%i.hdf5", nn);
#if OPERATOR==_MEHRSTELLEN_
            amr_op.vcycle(Phi, RhsCorr, Res, nn);
#else 
            amr_op.vcycle(Phi, Rhs, Res, nn);
#endif
            cout << scientific << "Residual: Max = " << absMax(Res);
            cout << "\t\tIntegral(Res) = " << integrate(Res, cdx) << endl;
        }
        amr_op.write(Res, "AMR_Res.%i.hdf5", numIter);
        amr_op.write(Phi, "AMR_Phi.%i.hdf5", numIter);

        for (int level = 0; level < numLevels; ++level)
        {
            auto& rhs = *(Rhs[level]);
            auto& sln = *(Sln[level]);
            auto& phi = *(Phi[level]);
            auto& err = *(Err[level]);
            auto iter = phi.dataIterator();
            for (iter.reset(); iter.ok(); ++iter)
            {
                BD sln_i = sln[iter];
                BD phi_i = phi[iter];
                BD err_i = err[iter];
                BD rhs_i = rhs[iter];
                phi_i.copyTo(err_i);
                err_i -= sln_i;
            }
        }
        amr_op.write(Err, "AMR_Error.hdf5");
        Real innerError = absMax(Err, 0, true, 1);
        Real bdryError = absMax(Err, 0, false, 1);

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
                Proto::Box domain = Proto::Box(Proto::Point(domainSize/s, domainSize));
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
                Proto::Box domain = Proto::Box(Proto::Point(domainSize, domainSize/s));
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
                std::cout << "\t\tDomain of level " << ii << ": " << domain << std::endl;
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
            AMRFAS<MehrstellenCorrectionOp<DATA>,DATA> correctOp(Layouts, dx*AMR_REFRATIO, 1);
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

            AMRFAS<OP> amr_op(Layouts, dx*AMR_REFRATIO , 1);

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
        LevelData<DATA> SlnC(coarseLayout, OP::numcomps(), Proto::Point::Zeros());
        LevelData<DATA> ResC(coarseLayout, OP::numcomps(), IntVect::Zero);
        LevelData<DATA> Phi(fineLayout, OP::numcomps(), IntVect::Unit);
        LevelData<DATA> Res(fineLayout, OP::numcomps(), IntVect::Zero);
        LevelData<DATA> Error(coarseLayout, OP::numcomps(), Proto::Point::Zeros());
#if OPERATOR==_MEHRSTELLEN_
        LevelData<DATA> RC(coarseLayout, OP::numcomps(), Proto::Point::Ones());
        LevelData<DATA> RCCorr(coarseLayout, OP::numcomps(), Proto::Point::Zeros());
        LevelData<DATA> R(fineLayout, OP::numcomps(), Proto::Point::Ones());
#else
        LevelData<DATA> RC(coarseLayout, OP::numcomps(), Proto::Point::Zeros());
        LevelData<DATA> R(fineLayout, OP::numcomps(), Proto::Point::Zeros());
#endif
        auto fiter = fineLayout.dataIterator();
        for (fiter.begin(); fiter.ok(); ++fiter)
        {
            BD res = Res[fiter];
            BD phi = Phi[fiter];
            BD rhs = R[fiter];
            res.setVal(0);
            phi.setVal(0);
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
            BD phiC = PhiC[citer];
            BD rhsC = RC[citer];
            BD errC = Error[citer];
            BD slnC = SlnC[citer];
            resC.setVal(0);
            phiC.setVal(0);
            errC.setVal(0);
            Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, Proto::Var<Real>& a_rhs)
            {
                Real x0 = a_p[0]*cdx;
                Real x1 = x0 + cdx;
                a_rhs(0) = (sin(x1) - sin(x0))/cdx;
            }, rhsC);
            Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, Proto::Var<Real>& a_sln)
            {
                Real x0 = a_p[0]*cdx;
                Real x1 = x0 + cdx;
                a_sln(0) = -(sin(x1) - sin(x0))/cdx;
            }, slnC);
#if OPERATOR==_MEHRSTELLEN_
            BD rhsCorr = RCCorr[citer];
            auto S = 1.0*Proto::Shift::Zeros() + 1.0/12.0*Proto::Stencil<Real>::Laplacian();
            rhsCorr |= S(rhsC);
#endif
        }

        writeLevel(RC, "MG_Rhs.hdf5");
        std::cout << "Integral of Rhs: " << integrate(RC, cdx) << std::endl;
#if OPERATOR==_MEHRSTELLEN_
        writeLevel(RCCorr, "MG_RhsCorr.hdf5");
        std::cout << "Integral of Corrected Rhs: " << integrate(RCCorr, cdx) << std::endl;
#endif

        if (doAMR)
        {
            std::cout << "AMR test is not fully functional" << std::endl;
            Multigrid<OP, DATA> amr_mg(fineLayout, dx, AMR_REFRATIO/MG_REFRATIO, true, 1);
            amr_mg.vcycle(Phi,PhiC,R);
        } else {
            Multigrid<OP, DATA> mg(coarseLayout, cdx, std::log2(domainSize), false);

            double resnorm = 0.0;
            OP op;
            op.define(coarseLayout,cdx);
            for (int ii = 0; ii < numIter; ii++)
            {
                writeLevel(PhiC, "MG_Phi.%i.hdf5",ii);
#if OPERATOR==_MEHRSTELLEN_
                mg.vcycle(PhiC, RCCorr);
                resnorm = op.residual(ResC,PhiC,RCCorr);
#else
                mg.vcycle(PhiC,RC); 
                resnorm = op.residual(ResC,PhiC,RC);
#endif
                std::cout << scientific << "iteration number = "
                    << ii << ", Residual norm: " << resnorm << std::endl;
            }
            auto iter = PhiC.dataIterator();
            for (iter.begin(); iter.ok(); ++iter)
            {
                OP::patch phi = PhiC[iter];
                OP::patch sln = SlnC[iter];
                OP::patch err = Error[iter];
                phi.copyTo(err);
                err -= sln;
            }
            writeLevel(Error, "MG_Error.hdf5");
            Real error = absMax(Error);
            std::cout << scientific << "Error: " << error << std::endl;
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
        LevelData<DATA> LPhiCInner(coarseFineLayout, OP::numcomps(), Proto::Point::Zeros());
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
                    Real y0 = a_p[1]*cdx;
                    Real y1 = y0 + cdx;
                    a_data(0) = (sin(x1) - sin(x0))/cdx;
                    a_data(0) *= (sin(y1) - sin(y0))/cdx;
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
                    Real y0 = a_p[1]*dx;
                    Real y1 = y0 + dx;
                    a_data(0) = (sin(x1) - sin(x0))/dx;
                    a_data(0) *= (sin(y1) - sin(y0))/dx;
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
                    //a_phiC(0) = 0;
                }
            }, flux, LphiC);
        }

        LPhiC.copyTo(LPhiCInner);
        
        writeLevel(PhiC, "RefluxTest_PhiC.hdf5");
        writeLevel(Phi, "RefluxTest_Phi.hdf5");
        writeLevel(Flux, "RefluxTest_Soln.hdf5");

        Real intFlux = integrate(Flux, cdx);
        Real intCoarse = integrate(LPhiC, cdx);
        Real intFine = integrate(LPhi, dx);
        Real intCoarseInner = integrate(LPhiCInner, cdx);

        std::cout << "Integral of L(Phi) on fine grid: " << intFine << std::endl;
        //std::cout << "Integral of L(PhiC) on coarse grid: " << intCoarse << std::endl;
        std::cout << "Integral of L(PhiC) on refined grid: " << intCoarseInner << std::endl;
        std::cout << "Int(L(Phi)) - Int(L(PhiC)): " << intFine - intCoarseInner << std::endl;
        std::cout << "Integral of composite register: " << intFlux << std::endl;
        std::cout << "Error: " << std::abs(intFine - intCoarseInner - intFlux) << std::endl;
    }
    else if (TEST == 6) {
        std::cout << "Testing Coarse Residual" << std::endl;
        Real L = 2.0*M_PI;

        Real error[numIter];
        Real boundError[numIter];
        for (int nn = 0; nn < numIter; nn++)
        {
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
            LevelData<DATA> Phi(fineLayout, OP::numcomps(), OP::ghost());
            LevelData<DATA> RhsC(coarseLayout, OP::numcomps(), Proto::Point::Zeros());
            LevelData<DATA> Rhs(fineLayout, OP::numcomps(), Proto::Point::Zeros());
            LevelData<DATA> ResC(coarseLayout, OP::numcomps(), Proto::Point::Zeros());
            auto citer = coarseLayout.dataIterator();
            for (citer.begin(); citer.ok(); ++citer)
            {
                OP::patch phiC = PhiC[citer];
                OP::patch rhsC = RhsC[citer];
                OP::patch resC = ResC[citer];
                Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, OP::var& a_data)
                        {
                        Real x0 = a_p[0]*cdx;
                        Real x1 = x0 + cdx;
                        a_data(0) = -(sin(x1) - sin(x0))/cdx;
                        }, phiC);
                phiC.copyTo(rhsC);
                rhsC *= -1;
                resC.setVal(0);
            }
            auto fiter = fineLayout.dataIterator();
            for (fiter.begin(); fiter.ok(); ++fiter)
            {
                OP::patch phi = Phi[fiter];
                OP::patch rhs = Rhs[fiter];
                Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, OP::var& a_data)
                        {
                        Real x0 = a_p[0]*dx;
                        Real x1 = x0 + dx;
                        a_data(0) = -(sin(x1) - sin(x0))/dx;
                        }, phi);
                phi.copyTo(rhs);
                rhs *= -1;
            }

            LevelFluxRegister Register;
            Register.define(fineLayout, coarseLayout, fineLayout.physDomain(), AMR_REFRATIO, OP::numcomps());

            OP op;
            op.define(fineLayout, dx, true);
            op.coarseResidual(ResC, RhsC, PhiC, Rhs, Phi, Register);

            std::cout << "Integral of Coarse Residual: " << integrate(ResC, cdx) << std::endl;
            writeLevel(PhiC, "Test_CoarseResidual_PhiC.hdf5");
            writeLevel(RhsC, "Test_CoarseResidual_RhsC.hdf5");
            writeLevel(Phi, "Test_CoarseResidual_Phi.hdf5");
            writeLevel(Rhs, "Test_CoarseResidual_Rhs.hdf5");
            writeLevel(ResC, "Test_CoarseResidual_ResC.hdf5");
            
            error[nn] = absMax(ResC);
            std::cout << "Error: " << error[nn] << std::endl;
            domainSize *= 2;
        }
        for (int ii = 1; ii < numIter; ii++)
        {
            std::cout << "Rate: " << log2(error[ii-1]/error[ii]) << std::endl;
        }
    }
#if CH_MPI
    CH_TIMER_REPORT();
    MPI_Finalize();
#endif
}
