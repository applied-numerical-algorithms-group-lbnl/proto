#include "ProtoAMR.H"
#include "InputParser.H"
#include "BoxOp_Laplace.H"
#include "TestFunc.H"

using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    using Proto::pout;
    HDF5Handler h5;

    int domainSize = 64;
    int boxSize = 16;
    int numIter = 3;
    int numLevels = 2;
    double physDomainSize = 1;
    int refRatio = 4;
    int testNum = 0;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    InputArgs args;
    args.add("domainSize", domainSize);
    args.add("boxSize",    boxSize);
    args.add("numIter",    numIter);
    args.add("testNum",    testNum);
    args.add("numLevels",  numLevels);
    args.add("periodic_x", periodicity[0]);
    args.add("periodic_y", periodicity[1]);
#if DIM > 2
    args.add("periodic_z", periodicity[2]);
#endif
    args.parse(argc, argv);
    args.print();
    
    if (testNum == 0 && procID() == 0)
    {
        std::cout << "TEST 0: MEMTYPE == HOST" << std::endl;
    } else if (testNum == 1 && procID() == 0)
    {
        std::cout << "TEST 1: MEMTYPE == DEVICE" << std::endl;
    } else if (procID() == 0)
    {
        std::cout << "UNKNOWN TEST NUMBER: " << testNum << std::endl;
        return;
    }
    
    typedef BoxOp_Laplace<double, HOST> OP; 

    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = physDomainSize / domainSize;
        
        std::vector<DisjointBoxLayout> layouts;
        layouts.resize(numLevels);
                
        Point boxSizeV = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain domain(domainBox, periodicity);
        layouts[0].define(domain, boxSizeV);

        Box refinedRegion = domainBox;
        for (int lvl = 1; lvl < numLevels; lvl++)
        {
            Point prevSize = refinedRegion.high() - refinedRegion.low() + Point::Ones();
            refinedRegion = refinedRegion.grow(-prevSize / 4).refine(refRatio);
            
            Box refinedRegionPatches = refinedRegion.coarsen(boxSizeV);
            std::vector<Point> fineLayoutPatches;
            bool skip = true;
            for (auto iter = refinedRegionPatches.begin(); iter.ok(); ++iter)
            {
                if (skip)
                {
                    skip = false;
                    continue;
                }
                fineLayoutPatches.push_back(*iter);
            }
            ProblemDomain fineDomain = layouts[lvl-1].domain().refine(Point::Ones(refRatio));
            layouts[lvl].define(fineDomain, fineLayoutPatches, boxSizeV);
        }
        std::vector<Point> refRatios;
        refRatios.push_back(Point::Ones(refRatio));
        AMRGrid grid(layouts, refRatios, numLevels);
        
        if (testNum == 0)
        {
            AMROp<BoxOp_Laplace, double, HOST> op(grid, dx);
            AMRData<double, OP::numState(), HOST>      Phi(grid,      Point::Ones(3));
            AMRData<double, OP::numState(), HOST>      LPhi(grid,     Point::Zeros());
            AMRData<double, OP::numState(), HOST>      LPhiSln(grid,  Point::Zeros());
            AMRData<double, OP::numState(), HOST>      LPhiErr(grid,  Point::Zeros());
            LevelBoxData<double, OP::numState(), HOST> LPhi0(grid[0], Point::Zeros());
            LevelBoxData<double, OP::numState(), HOST> LPhiErr0(grid[0], Point::Zeros());
            LevelBoxData<double, OP::numState(), HOST> LPhi1(grid[1], Point::Zeros());
            LevelBoxData<double, OP::numState(), HOST> LPhiErr1(grid[1], Point::Zeros());

            Phi.initialize(dx, sinProd_avg);
            LPhiSln.initialize(dx, LsinProd_avg);

            LPhi.setToZero();
            LPhi0.setToZero();
            LPhi1.setToZero();
            LPhiErr.setToZero();
            LPhiErr0.setToZero(); 
            LPhiErr1.setToZero(); 

            Phi.averageDown();
            auto INTERP = InterpStencil<double>::Build(OP::order()+1, refRatio);
            interpBoundaries(Phi[0], Phi[1], INTERP);
            Phi.exchange();
            op(LPhi, Phi);
            //op.levelApply(LPhi0, Phi, 0);
            //op.levelApply(LPhi1, Phi, 1);
            LPhi.averageDown();
            
            LPhi.copyTo(LPhiErr);
            LPhiErr.increment(LPhiSln, -1);
            
            err[nn] = LPhiErr.absMax();
                    
            if (procID() == 0)
            {
                std::cout << "Error: " << err[nn] << std::endl;
            }
        } else if (testNum == 1)
        {
            AMROp<BoxOp_Laplace, double, DEVICE> op(grid, dx);
            AMRData<double, OP::numState(), DEVICE>      Phi(grid,      Point::Ones(3));
            AMRData<double, OP::numState(), DEVICE>      LPhi(grid,     Point::Zeros());
            AMRData<double, OP::numState(), DEVICE>      LPhiSln(grid,  Point::Zeros());
            AMRData<double, OP::numState(), DEVICE>      LPhiErr(grid,  Point::Zeros());
            LevelBoxData<double, OP::numState(), DEVICE> LPhi0(grid[0], Point::Zeros());
            LevelBoxData<double, OP::numState(), DEVICE> LPhiErr0(grid[0], Point::Zeros());
            LevelBoxData<double, OP::numState(), DEVICE> LPhi1(grid[1], Point::Zeros());
            LevelBoxData<double, OP::numState(), DEVICE> LPhiErr1(grid[1], Point::Zeros());

            Phi.initialize(dx, sinProd_avg);
            LPhiSln.initialize(dx, LsinProd_avg);

            LPhi.setToZero();
            LPhi0.setToZero();
            LPhi1.setToZero();
            LPhiErr.setToZero();
            LPhiErr0.setToZero(); 
            LPhiErr1.setToZero(); 

            Phi.averageDown();
            auto INTERP = InterpStencil<double>::Build(OP::order()+1, refRatio);
            interpBoundaries(Phi[0], Phi[1], INTERP);
            Phi.exchange();
            op(LPhi, Phi);
            //op.levelApply(LPhi0, Phi, 0);
            //op.levelApply(LPhi1, Phi, 1);
            LPhi.averageDown();
            
            LPhi.copyTo(LPhiErr);
            LPhiErr.increment(LPhiSln, -1);
            
            err[nn] = LPhiErr.absMax();
                    
            if (procID() == 0)
            {
                std::cout << "Error: " << err[nn] << std::endl;
            }
        }
        domainSize *= 2;
    }
        
    for (int ii = 1; ii < numIter; ii++)
    {
        std::cout << "Convergence Rate: " << log(err[ii-1] / err[ii]) / log(2.0) << std::endl;
    }

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

