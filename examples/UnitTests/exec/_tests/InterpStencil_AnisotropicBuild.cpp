#include "Proto.H"
#include "InputParser.H"
#include "func.H"

using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif

    // SETUP
    using Proto::pout;
    HDF5Handler h5;
    InputArgs args;
    args.parse();
    args.print();

    int domainSize = 64;
    int boxSize = 16;
    int numIter = 3;
    int ghostSize = 1;
    int refRatio = 2;
    int order = 5;
    double k = 1;
    double physDomainSize = 1;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    args.set("domainSize", &domainSize);
    args.set("boxSize",    &boxSize);
    args.set("numIter",    &numIter);
    args.set("refRatio",   &refRatio);
    args.set("order",      &order);
    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);

    Point refineX(refRatio,1,1,1,1,1);
    Point refineY(1,refRatio,1,1,1,1);
    Point refineXY(refRatio,refRatio,1,1,1,1);
    
    double err_x[numIter];
    double err_y[numIter];
    double err_xy[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = 1.0/domainSize;
        std::array<double, DIM> fdx_x;
        std::array<double, DIM> fdx_y;
        std::array<double, DIM> fdx_xy;
        for (int dir = 0; dir < DIM; dir++)
        {
            fdx_x[dir] = dx / refineX[dir];
            fdx_y[dir] = dx / refineY[dir];
            fdx_xy[dir] = dx / refineXY[dir];
        }

        auto interpX  = InterpStencil<double>::Build(order, refineX);
        auto interpY  = InterpStencil<double>::Build(order, refineY);
        auto interpXY = InterpStencil<double>::Build(order, refineXY);

        interpX.print();

        Box domainBox  = Box::Cube(domainSize);
        Box rangeBoxX  = domainBox.refine(refineX);
        Box rangeBoxY  = domainBox.refine(refineY);
        Box rangeBoxXY = domainBox.refine(refineXY);
        
        BoxData<double> input(domainBox.grow(interpXY.ghost()));
        
        BoxData<double> outputX(rangeBoxX);
        BoxData<double> outputY(rangeBoxY);
        BoxData<double> outputXY(rangeBoxXY);
        
        BoxData<double> solutionX(rangeBoxX);
        BoxData<double> solutionY(rangeBoxY);
        BoxData<double> solutionXY(rangeBoxXY);
        
        BoxData<double> errorX(rangeBoxX);
        BoxData<double> errorY(rangeBoxY);
        BoxData<double> errorXY(rangeBoxXY);

        forallInPlace_p(f_wave_avg, input, dx,  k);
        forallInPlace_p(f_wave_avg_aniso, solutionX,  fdx_x,  k);
        forallInPlace_p(f_wave_avg_aniso, solutionY,  fdx_y,  k);
        forallInPlace_p(f_wave_avg_aniso, solutionXY, fdx_xy, k);

        outputX.setVal(0);
        outputY.setVal(0);
        outputXY.setVal(0);
        
        h5.writePatch(dx, input, "INPUT_N%i", nn);
        outputX  |= interpX(input);
        outputY  |= interpY(input);
        outputXY |= interpXY(input);
        h5.writePatch(fdx_x,  outputX,  "OUTPUT_X_N%i", nn);
        h5.writePatch(fdx_y,  outputY,  "OUTPUT_Y_N%i", nn);
        h5.writePatch(fdx_xy, outputXY, "OUTPUT_XY_N%i", nn);
        h5.writePatch(fdx_x,  solutionX,  "SOLN_X_N%i", nn);
        h5.writePatch(fdx_y,  solutionY,  "SOLN_Y_N%i", nn);
        h5.writePatch(fdx_xy, solutionXY, "SOLN_XY_N%i", nn);
        
        outputX.copyTo(errorX);
        outputY.copyTo(errorY);
        outputXY.copyTo(errorXY);
        errorX  -= solutionX;
        errorY  -= solutionY;
        errorXY -= solutionXY;
        h5.writePatch(fdx_x,  errorX,  "ERROR_X_N%i", nn);
        h5.writePatch(fdx_y,  errorY,  "ERROR_Y_N%i", nn);
        h5.writePatch(fdx_xy, errorXY, "ERROR_XY_N%i", nn);
        err_x[nn] = errorX.absMax();
        err_y[nn] = errorY.absMax();
        err_xy[nn] = errorXY.absMax();

        pout() << "Error (X): " << err_x[nn] << std::endl;
        pout() << "Error (Y): " << err_y[nn] << std::endl;
        pout() << "Error (XY): " << err_xy[nn] << std::endl;
        domainSize *= 2;
    }
        
    for (int ii = 1; ii < numIter; ii++)
    {
        pout() << "Convergence Rate (X): "  << log(err_x[ii-1] / err_x[ii]) / log(2.0) << std::endl;
        pout() << "Convergence Rate (Y): "  << log(err_y[ii-1] / err_y[ii]) / log(2.0) << std::endl;
        pout() << "Convergence Rate (XY): " << log(err_xy[ii-1] / err_xy[ii]) / log(2.0) << std::endl;
    }

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

