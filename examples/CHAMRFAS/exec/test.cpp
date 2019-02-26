#include "LevelData.H"
#include "Multigrid.H"
#include "Proto.H"
#include "TestOp.H" //definition of OP
#include "AMRFAS.H"

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
        cout << "\t\tTest 0: Interp Testing" << endl;
        cout << "\t\tTest 1: AMR Multigrid" << endl;
        cout << "\t\tTest 2: AMR Operator" << endl;
        cout << "\t\tTest 3: Multigrid" << endl;
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

    typedef Proto::BoxData<Real, NUMCOMPS> BD;
    typedef TestOp<FArrayBox> OP;
    typedef FArrayBox DATA;

     
    //====================================================================
    // Misc Testing
    if (TEST == 0)
    {
        //Real L = 2.0*M_PI;
        Real L = domainSize;
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

            OP op(fineLayout, dx);

            LevelData<DATA> PhiF(fineLayout, 1, Proto::Point::Ones());
            LevelData<DATA> FluxX(fineLayout, 1, Proto::Point::Ones());
            LevelData<DATA> TrueFluxX(fineLayout, 1, Proto::Point::Ones());
            LevelData<DATA> FluxY(fineLayout, 1, Proto::Point::Ones());
            LevelData<DATA> AvgFluxX(coarseTemp, 1, Proto::Point::Ones());
            LevelData<DATA> AvgFluxY(coarseTemp, 1, Proto::Point::Ones());
            LevelData<DATA> LPhiF(fineLayout, 1);
            for (fiter.begin(); fiter.ok(); ++fiter)
            {
                BD data = PhiF[fiter];
                BD trueFlux = TrueFluxX[fiter];
                data.setVal(0);
                Proto::Box patchBox = fineLayout[fiter];
                Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, V& a_data)
                        {
                        Real x0 = dx*a_p[0];
                        Real x1 = x0 + dx;
                        a_data(0) = (sin(x1) - sin(x0))/dx;
                        //a_data(0) = 4;
                        }, patchBox, data);
                Proto::forallInPlace_p([=] PROTO_LAMBDA (Proto::Point& a_p, V& a_data)
                        {
                            a_data(0) = -sin(dx*a_p[0]);
                        },trueFlux);
            }
            PhiF.exchange();
            for (fiter.reset(); fiter.ok(); ++fiter)
            {
                BD phi = PhiF[fiter];
                BD lphi = LPhiF[fiter];
                BD flux_x = FluxX[fiter];
                BD flux_y = FluxY[fiter];
                op(lphi, phi);
            }
            op.flux(FluxX, PhiF, 0);
            op.flux(FluxY, PhiF, 1);
            FluxX.exchange();
            FluxY.exchange();
            
            auto avgx = 0.5*Proto::Shift::Zeros() + 0.5*Proto::Shift::Basis(0);
            avgx.srcRatio() = Proto::Point::Ones(2);
            auto avgy = 0.5*Proto::Shift::Zeros() + 0.5*Proto::Shift::Basis(1);
            avgy.srcRatio() = Proto::Point::Ones(2);

            for (fiter.reset(); fiter.ok(); ++fiter)
            {
                BD fluxX = FluxX[fiter];
                BD avgFluxX = AvgFluxX[fiter];
                avgFluxX.setVal(0);
                avgFluxX |= avgy(fluxX,coarseTemp[fiter], 1.0/cdx);
                
                BD fluxY = FluxY[fiter];
                BD avgFluxY = AvgFluxY[fiter];
                avgFluxY.setVal(0);
                avgFluxY |= avgx(fluxY, coarseTemp[fiter], 1.0/cdx);
            }

            LevelData<DATA> FluxRegister(coarseLayout, 1, Proto::Point::Zeros()); //surrogate for RC
            LevelData<DATA> PhiC(coarseLayout, 1, Proto::Point::Ones());
            for (citer.begin(); citer.ok(); ++citer)
            {
                BD fr = FluxRegister[citer];
                BD phic = PhiC[citer];
                fr.setVal(0);
                phic.setVal(0);
            }
            LevelFluxRegister LFR;
            LFR.define(fineLayout, coarseLayout, fineLayout.physDomain(), AMR_REFRATIO, 1);

            writeLevel(PhiF, "PhiF.test.hdf5");
            writeLevel(LPhiF, "LPhiF.test.hdf5");
            writeLevel(PhiC, "PhiC.test.hdf5");
            writeLevel(FluxX, "FluxX.test.hdf5");
            writeLevel(TrueFluxX, "TrueFluxX.test.hdf5");
            writeLevel(FluxY, "FluxY.test.hdf5");
            writeLevel(AvgFluxX, "AvgFluxX.test.hdf5");
            writeLevel(AvgFluxY, "AvgFluxY.test.hdf5");

            op.reflux(FluxRegister, PhiC, PhiF, LFR);

            for (citer.reset(); citer.ok(); ++citer)
            {
                BD data = FluxRegister[citer];
                Proto::Box patchBox = coarseLayout[citer];
                Proto::Box intersect = patchBox & coarseFineDomain;
                Proto::forallInPlace([=] PROTO_LAMBDA (V& a_data)
                        {
                        a_data(0) = 0;
                        }, intersect, data);
            }
            writeLevel(FluxRegister, "FluxRegister.test.hdf5");

            Real lphi = integrate(LPhiF, dx);
            Proto::Box halo = coarseFineDomain.grow(1);
            Real flux = 0.0;
            for (citer.reset(); citer.ok(); ++citer)
            {
                BD data = FluxRegister[citer];
                Proto::Box b = coarseLayout[citer];
                for (auto iter = b.begin(); iter != b.end(); ++iter)
                {
                    if (halo.contains(*iter) && (!coarseFineDomain.contains(*iter)))
                    {
                        flux += data(*iter)*cdx*cdx;
                    }
                }
            }

            //Real flux = integrate(FluxRegister, cdx);

            cout << "Integral of LPhiF: " << integrate(LPhiF, dx) << endl;
            cout << "Integral of FluxRegister: " << integrate(FluxRegister, cdx) << endl;
            cout << "Error: " << flux + lphi << endl;
            cout << "cdx/Error: " << cdx/(flux+lphi) << endl;
            error[n] = flux + lphi;
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
        cout << pow(0.5,DIM) << endl;
        int numLevels = 1;
        bool do_nesting = false;

        Real cdx = 2.0*M_PI/domainSize;
        std::vector<DisjointBoxLayout> Layouts;
        
        std::vector<std::shared_ptr<LevelData<DATA>>> Phi;
        std::vector<std::shared_ptr<LevelData<DATA>>> Src;
        std::vector<std::shared_ptr<LevelData<DATA>>> Rhs;
        std::vector<std::shared_ptr<LevelData<DATA>>> Res;
        std::vector<std::shared_ptr<LevelData<DATA>>> Const;

        Layouts.resize(numLevels);
        Phi.resize(numLevels);
        Src.resize(numLevels);
        Rhs.resize(numLevels);
        Res.resize(numLevels);
        Const.resize(numLevels);

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
                std::cout << "Domain of level " << ii << ": " << domain << std::endl;
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
                std::cout << "Domain of level " << ii << ": " << domain << std::endl;
                buildLayout(layout, domain, Proto::Point::Ones());
            }
            Phi[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Ones(1));
            Src[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Ones(1));
            Rhs[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
            Res[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Zeros());
            Const[ii] = std::make_shared<LevelData<DATA>>(layout, 1, Proto::Point::Ones());
            auto& phi = *(Phi[ii]);
            auto& src = *(Src[ii]);
            auto& rhs = *(Rhs[ii]);
            auto& res = *(Res[ii]);
            auto& con = *(Const[ii]);
            auto iter = phi.dataIterator();
            for (iter.reset(); iter.ok(); ++iter)
            {
                BD phi_i = phi[iter];
                BD rhs_i = rhs[iter];
                BD res_i = res[iter];
                BD con_i = con[iter];
                con_i.setVal(3);
                res_i.setVal(0);
                forallInPlace_p(
                        [=] PROTO_LAMBDA (Proto::Point& a_pt, Proto::Var<Real>& a_phi)
                        {
                        //Real x0 = a_pt[0]*dx;
                        //Real x1 = a_pt[0]*dx + dx;
                        //a_phi(0) = -0.5*(sin(x1) - sin(x0))/dx;
                        //a_phi(0) += 0.1*(sin(5*x1)/5 - sin(5*x0)/5)/dx;
                        a_phi(0) = 0.0; 
                        }, phi_i);
                forallInPlace_p(
                        [=] PROTO_LAMBDA (Proto::Point& a_pt, Proto::Var<Real>& a_rhs)
                        {
                        Real x0 = a_pt[0]*dx;
                        Real x1 = a_pt[0]*dx + dx;
                        a_rhs(0) = (sin(x1) - sin(x0))/dx;
                        }, rhs_i);
            }
            phi.exchange();
            phi.copyTo(src);
            dx /= AMR_REFRATIO;
        } //end initialization
        
        cout << "Error of Integral of Const: " << integrate(Const, cdx) - pow(2*M_PI, DIM)*3.0 << endl;

        AMRFAS<OP,DATA> amr_op(Layouts, dx*AMR_REFRATIO , numLevels-1, 1);
        amr_op.write(Src, "AMR_Src.hdf5");
        cout << "Integral of RHS: " << integrate(Rhs, cdx) << endl; 
        cout << "Integral of Res: " << integrate(Res, cdx) << endl; 
        cout << "Integral of Phi: " << integrate(Phi, cdx) << endl; 
        for (int nn = 0; nn < numIter; nn++)
        {
            amr_op.write(Res, "AMR_Res.%i.hdf5", nn);
            amr_op.write(Phi, "AMR_Phi.%i.hdf5", nn);
            amr_op.vcycle(Phi, Rhs, Res);
            cout << "Residual: Max = " << absMax(Res) << " Integral = " << integrate(Res, cdx) << endl;
        }
        amr_op.write(Res, "AMR_Res.%i.hdf5", numIter);
        amr_op.write(Phi, "AMR_Phi.%i.hdf5", numIter);

        amr_op.write(Rhs, "AMR_Rhs.hdf5");
        amr_op.write(Phi, "AMR_Phi.hdf5");
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
                rhsPatch += phiPatch;
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
        const int numLevels = 2;
        Real ei[numIter];
        Real eb[numIter];
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
#if CH_MPI 
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
        int numLevels = log(domainSize*1.0)/log(2.0)-1;
        Real dx = 2.0*M_PI/domainSize;
#if CH_MPI
        if (rank == 0)
        { 
#endif
            std::cout << "Running Multigrid:" << std::endl;
            std::cout << "\tDomain Size: " << domainSize << std::endl;
            std::cout << "\tMax Box Size: " << MAXBOXSIZE << std::endl;
            std::cout << "\tNumber of Multigrid Levels: " << numLevels << std::endl;
#if CH_MPI
        } 
#endif
        Box domainBox = Proto::Box::Cube(domainSize);
        DisjointBoxLayout layout;
        buildLayout(layout, domainBox);

        LevelData<DATA> U(layout, NUMCOMPS, IntVect::Unit);
        LevelData<DATA> F(layout, NUMCOMPS, IntVect::Zero);
        LevelData<DATA> R(layout, NUMCOMPS, IntVect::Zero);
        LevelData<DATA> S(layout, NUMCOMPS, IntVect::Zero);

        OP::initialCondition(U,dx);
        OP::forcing(F,dx);
        OP::solution(S,dx);

        Multigrid<OP, DATA> mg(layout, dx, numLevels-1, false);
        int numIter = 20;
        double resnorm = 0.0;
        //char fileName[100];
        //char fileNameU[100];
        TestOp<FArrayBox> op(layout,dx);
        //fileNum = 0;
        for (int ii = 0; ii < numIter; ii++)
        {
            mg.vcycle(U,F); 
            resnorm = op.residual(R,U,F);
#if CH_MPI
            if (rank == 0)
            {
#endif
                std::cout << scientific << "iteration number = " << ii << ", Residual norm: " << resnorm << std::endl;
#if CH_MPI
            }
#endif
            /*
            sprintf(fileName,"ResV.%i.hdf5",fileNum);
            sprintf(fileNameU,"ResU.%i.hdf5",fileNum);
            writeLevelname(&R,fileName);
            writeLevelname(&U,fileNameU);
            fileNum++;
            */
        }
        auto iter = layout.dataIterator();
        double umax = 0.0;
        double umin = 0.0;
        double error = 0.0;
        for (iter.begin(); iter.ok(); ++iter)
        {
            Proto::BoxData<double> s = S[iter()];
            Proto::BoxData<double> u = U[iter()];
            umax = max(umax,u.max());
            umin = min(umin,u.min());
            s -= u;
            error = max(s.absMax(),error);
        }
#if CH_MPI
        double max_error;
        MPI_Reduce(&error, &max_error, 1,  MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (rank == 0)
        {
            cout << "Error: " << max_error << endl;
#else
            cout << "Error: " << error << endl;
#endif
#if CH_MPI
        } 
#endif
    } // End Multigrid test
#if CH_MPI
    CH_TIMER_REPORT();
    MPI_Finalize();
#endif
}
