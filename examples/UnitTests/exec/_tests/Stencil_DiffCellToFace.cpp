#include "Proto.H"
#include "func.H"

using namespace Proto;

int main(int argc, char** argv)
{
    HDF5Handler h5;
    double k = 1;

    std::vector<Stencil<double>> DIFF;
    Point span = Point::Zeros();
    for (int dir = 0; dir < DIM; dir++)
    {
        auto D = Stencil<double>::DiffCellToFace(dir);
        DIFF.push_back(D);
        span += DIFF[dir].ghost();
    }

    int domainSize = 32;
    int numIter = 3;
    double err[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = 1.0 / domainSize;
        Box domain = Box::Cube(domainSize);

        BoxData<double> phi(domain.grow(span));
        forallInPlace_p(f_wave_avg, phi, dx, k);
        h5.writePatch(dx, phi, "PHI_N%i", nn);

        BoxData<double, DIM> dphi_soln(domain);
        forallInPlace_p(f_dwave_avg, dphi_soln, dx, k);
        h5.writePatch(dx, dphi_soln, "DPHI_SOLN_N%i", nn);
        
        BoxData<double, DIM> dphi(domain);
        dphi.setVal(7);

        for (int dir = 0; dir < DIM; dir++)
        {
            auto dphi_i = slice(dphi, dir);
            dphi_i |= DIFF[dir](phi, 1.0/dx);
        }

        BoxData<double, DIM> error(domain);
        dphi.copyTo(error);
        error -= dphi_soln;

        h5.writePatch(dx, dphi, "DPHI_N%i", nn);
        h5.writePatch(dx, error, "DPHI_ERROR_N%i", nn);

        err[nn] = error.absMax();
        std::cout << "Error: " << err[nn] << std::endl;
        domainSize *= 2;
    }

    for (int ii = 1; ii < numIter; ii++)
    {
        std::cout << "Convergence Rate: " << log(err[ii-1]/err[ii])/log(2.0) << std::endl;
    }
}
