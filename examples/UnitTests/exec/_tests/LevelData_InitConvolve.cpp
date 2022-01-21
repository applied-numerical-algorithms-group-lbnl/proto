#include "Proto.H"
#include "../../common/InputParser.H"
#include "func.H"

using namespace Proto;

int main(int argc, char** argv)
{
    HDF5Handler h5;

    int domainSize = 32;
    int boxSize = 8;
    int ghostSize = 0;
    int numIter = 3;
    std::array<bool, DIM> periodicity;

    InputArgs args;
    args.parse();
    args.set("domainSize", &domainSize);
    args.set("boxSize", &boxSize);
    args.set("numIter", &numIter);

    double k = 1;
    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = 1.0 / domainSize;

        Point boxSizeV = Point::Ones(boxSize);
        Box domainBox = Box::Cube(domainSize);
        ProblemDomain domain(domainBox, periodicity);
        DisjointBoxLayout layout(domain, boxSizeV);

        LevelBoxData<double> analytic(layout, Point::Ones(ghostSize));
        analytic.initialize(f_wave_avg, dx, k); 
        LevelBoxData<double> convolved(layout, Point::Ones(ghostSize));
        convolved.initConvolve(f_wave, dx, k); 
        LevelBoxData<double> error(layout, Point::Zeros());
        error.setToZero();

        for (auto iter = layout.begin(); iter.ok(); ++iter)
        {
            auto& analytic_i  = analytic[*iter];
            auto& convolved_i = convolved[*iter];
            auto& error_i = error[*iter];
            convolved_i.copyTo(error_i);
            error_i -= analytic_i;
        }
        err[nn] = error.absMax();
        std::cout << "Error: " << err[nn] << std::endl;
        domainSize *= 2;
    }
    for (int ii = 1; ii < numIter; ii++)
    {
        std::cout << "Convergence Rate: " << log(err[ii-1]/err[ii])/log(2.0) << std::endl;
    }
}
