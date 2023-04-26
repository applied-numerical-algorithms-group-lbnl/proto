#include <iostream>
#include <iomanip>
#include "Kokkos_Core.hpp"
#include "Proto.H"
#include "BoxOp_Advection.H"
#include "AMRRK4.H"
#include "InputParser.H"

using namespace Proto;

int TIME_STEP = 0;
int PLOT_NUM = 0;
PROTO_KERNEL_START
void
f_initializeF(
        const Point& a_pt,
        Var<double,NUMCOMPS>& a_U,
        const double& a_h,
        const double& a_time)
{
    for (int comp = 0; comp < NUMCOMPS; comp++)
    {
        a_U(comp) = 1.;
        for (int dir = 0; dir < DIM ; dir++)
        {
            a_U(comp)*=(-cos(2*M_PI*(a_pt[dir]*a_h - a_time) + 2*M_PI*a_h) 
                      + cos(2*M_PI*(a_pt[dir]*a_h - a_time)))/(2*M_PI*a_h);
        }
    }
}
PROTO_KERNEL_END(f_initializeF, f_initialize);

PROTO_KERNEL_START
void f_advectionExactF(
        Point& a_pt,
        Var<double,NUMCOMPS>& a_U,
        double a_h,
        double a_time)
{
    double r0 = .125;
    for (int comp = 0; comp < NUMCOMPS; comp++)
    {
        double rsq = 0.;
        for (int dir = 0; dir < DIM ; dir++)
        {
            double xcen = fmod(.5 + a_time,1.); 
            double xdir = a_pt[dir]*a_h + .5*a_h;
            double del; 
            if (xcen > xdir) del = min(abs(xdir - xcen),abs(xdir - xcen + 1.));
            if (xcen <= xdir) del = min(abs(xdir - xcen),abs(xdir - xcen - 1.));
            rsq += pow(del,2);
        }
        double r = sqrt(rsq);
        if (r > r0)
        {
            a_U(comp) = 0.;
        }
        else
        {
            a_U(comp) = pow(cos(M_PI*r/r0/2),6);
        }
    }
}
PROTO_KERNEL_END(f_advectionExactF, f_advectionExact);

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc,&argv);
#endif
  Kokkos::initialize(argc, argv);
  {
#ifdef PR_HDF5
    HDF5Handler h5;
#endif
    using Proto::pout; 
    typedef BoxOp_Advection<double> OP;

    int domainSize =32;
    int boxSize = 32;
    int refRatio = 4;
    int maxTimesteps = 4096;
    int numLevels = 3;
    int regridBufferSize = 3;
    int regridInterval = 1;
    double maxTime = 1.0;
    int outputInterval = 16;
    double t0 = 0.;
    int maxRefs = 2;
    //std::array<bool, DIM> periodicity;
    Array<bool, DIM> periodicity;
    periodicity.fill(true);

    InputArgs args;
    args.add("domainSize",      domainSize);
    args.add("boxSize",         boxSize);
    args.add("refRatio",        refRatio);
    args.add("maxTimesteps",    maxTimesteps);
    args.add("maxTime",         maxTime);
    args.add("outputInterval",  outputInterval);
    args.add("periodic_x",      periodicity[0]);
    args.add("periodic_y",      periodicity[1]);
    args.add("numLevels",       numLevels);
    args.add("regridBufferSize",regridBufferSize);
    args.add("regridInterval"  ,regridInterval);
    args.add("t0",              t0);
    args.add("maxRefs",         maxRefs);
    args.parse(argc, argv);
    std::vector<double> errorRefs;
    double consError;
    for (int numRef = 0;numRef < maxRefs;numRef++)
      {
        args.print();

        if (boxSize > domainSize)
          {
            pout() << "Input error: boxSize > domainSize. Forcing boxSize == domainSize." << std::endl;
            boxSize = domainSize;
          }
        PR_TIMER_SETFILE(to_string(domainSize) 
                         + ".DIM=" + to_string(DIM) + ".numProc=" + to_string(numProc())
                         + "AMRAdvection.time.table");

        ofstream sizeout(to_string(domainSize) 
                         + ".DIM=" + to_string(DIM) + ".numProc=" + to_string(numProc())
                         + "HWM.curve");
        Point boxSizeVect = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain domain(domainBox, periodicity);
        DisjointBoxLayout layout(domain, boxSizeVect);

        double physDomainSize = 1.0;
        //std::array<double, DIM> dx;
        Array<double, DIM> dx;
        dx.fill(physDomainSize / domainSize);
        double dt = 0.5 / domainSize;

        pout() << setw(50) << setfill('=') << "=" << std::endl;
        pout() << "Coarsest Level Parameters" << std::endl;
        pout() << setw(50) << setfill('-') << "-" << std::endl;
        pout() << "dx: " << dx[0] << " | dt: " << dt << std::endl;
        pout() << "domain: " << domainBox << std::endl;
        pout() << setw(50) << setfill('-') << "-" << std::endl << std::endl;

        std::vector<Point> refRatios(numLevels-1, Point::Ones(refRatio));
        AMRGrid grid(layout, refRatios, numLevels);

        double dxLevel = dx[0];
        Point bufferSize = Point::Ones(regridBufferSize);
        if (numLevels > 1)
          {
            for (int lvl = 0; lvl < numLevels-1; lvl++)
              {
                LevelBoxData<double, NUMCOMPS> initData(grid[lvl], Point::Zeros());
                Operator::initConvolve(initData, f_advectionExact, dxLevel, t0);
                LevelTagData tags(grid[lvl], bufferSize);
                LevelTagData kokkos_tags(grid[lvl], bufferSize);
                for (auto iter = grid[lvl].begin(); iter.ok(); ++iter)
                  {
                    //OP::generateTags(tags[*iter], initData[*iter]);
                  }
                for (auto iter = grid[lvl].begin(); iter.ok(); ++iter)
                  {
                    OP::kokkos_generateTags(kokkos_tags[*iter], initData[*iter]);
                  }
                AMRGrid::buffer(tags, bufferSize);
                grid.regrid(tags, lvl);
                dxLevel /= refRatio;
              }
        
            for (int lvl = 2; lvl < numLevels; lvl++)
              {
                grid.enforceNesting2(lvl);
              }
          }
        AMRData<double, NUMCOMPS> U(grid, OP::ghost());
        Operator::initConvolve(U, dx[0], f_advectionExact, t0);
        U.averageDown();
        AMRRK4<OP, double, NUMCOMPS> advectionOp(U, dx,regridInterval,regridBufferSize);
        double sum0 = U[0].integrate(dx[0]);
        pout() << "Initial level 0 conservation sum: " << sum0 << std::endl;

        double time = t0;
#ifdef PR_HDF5
        h5.setTime(time);
        h5.setTimestep(dt);
        h5.writeAMRData(dx, U, "U_"+to_string(domainSize)+"_N0");
#endif
        {      
          struct rusage usage;
          getrusage(RUSAGE_SELF, &usage);
          if (Proto::procID() == 0)
            {
              cout << "memory HWM before evolution: " << usage.ru_maxrss << endl;
              int step = 0;
              sizeout << step << " " << usage.ru_maxrss << endl;
            }
          // pout() << "memory HWM before evolution: " << usage.ru_maxrss << endl;
        }
        for (int k = 0; ((k < maxTimesteps) && (time < maxTime)); k++)
          {
            TIME_STEP = k;
            PLOT_NUM = 0;
         
            advectionOp.advance(dt);

            // Memory HWM.
        
            struct rusage usage;
            getrusage(RUSAGE_SELF, &usage);
            if (Proto::procID() == 0)
              {
                sizeout << k << " " << usage.ru_maxrss << endl;
              }      
            time += dt;
#ifdef PR_HDF5
            h5.setTime(time);
#endif
            if ((k+1) % outputInterval == 0)
              {
                int numBoxes = 0;
                for (int ll = 0; ll < U.numLevels();ll++)
                  {
                    numBoxes += U[ll].size();
                  }
                if ((Proto::procID() == 0)&&(((k+1) % outputInterval) == 0 ))
                  {
                    std::cout << "n: " << k+1 << " |  time: " << time << 
                      "  |  number of Boxes: " << numBoxes << std::endl;
                  }
#ifdef PR_HDF5             
                h5.writeAMRData(dx, U, "U_"+to_string(domainSize)+"_N"+to_string(k+1));
#endif
              }        
          }
        {
          struct rusage usage;
          getrusage(RUSAGE_SELF, &usage);
          if (Proto::procID() == 0)
            cout << "memory HWM after evolution: " << usage.ru_maxrss << endl;
          // pout() << "memory HWM after evolution: " << usage.ru_maxrss << endl;
        }
        U.averageDown();
        
        AMRData<double, NUMCOMPS> USoln(U.grid(), Point::Zeros());
        Operator::initConvolve(USoln, dx[0], f_advectionExact, time);

        AMRData<double, NUMCOMPS> UError(U.grid(), Point::Zeros());
        U.copyTo(UError);
        UError.increment(USoln, -1);
    
#ifdef PR_HDF5
        h5.writeAMRData(dx, UError, "UError");
#endif
        double l1Soln = USoln.integrateAbs(dx);
        double error = UError.integrateAbs(dx)/l1Soln;
        errorRefs.push_back(error);
        if (Proto::procID() == 0)
          {
            std::cout << "error: " << error << std::endl;
            std::cout << " time: " << time << std::endl;
          }
        double sumFinal = U[0].integrate(dx[0]);
        consError = (sumFinal - sum0)/l1Soln;
        if (Proto::procID() == 0)
          {
            //cout << "Final level 0 conservation sum: " << sumFinal << std::endl;   
            cout << "Relative conservation error: " << consError << std::endl;
            
          }
        domainSize *=2;
        dt *= .5;
        maxTimesteps *= 2;
      }
    if (Proto::procID() == 0)
      std::cout <<
        "Relative conservation error at finest level (must be < 1.e-14 to pass): "
                << fabs(consError) << std::endl;
    if (maxRefs > 1)
      {
        double rate = log(errorRefs[maxRefs-2]/errorRefs[maxRefs-1])/log(2.0);
        if (Proto::procID() == 0)
          {
            std::cout << "L1 Richardson error exponent (must be >= 3.9 to pass): "
                      << rate << std::endl;         
          }
      }
    PR_TIMER_REPORT();
    }
    Kokkos::finalize();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
