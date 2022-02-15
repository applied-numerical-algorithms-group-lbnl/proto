#include "Proto.H"
#include "func.H"

using namespace Proto;

int main(int argc, char** argv)
{
    HDF5Handler h5;
    double k = 1;

    std::vector<Stencil<double>> C2F;
    Point span = Point::Zeros();
    for (int dir = 0; dir < DIM; dir++)
    {
        auto D = Stencil<double>::CellToFace(dir, Side::Lo, 5);
        C2F.push_back(D);
        span += C2F[dir].ghost();
        std::cout << "Default Span: " << D.span() << std::endl;
    }
    std::vector<Stencil<double>> C2F_L;
    Point span_L = Point::Zeros();
    for (int dir = 0; dir < DIM; dir++)
    {
        auto D = Stencil<double>::CellToFaceL(dir, Side::Lo, 5);
        C2F_L.push_back(D);
        span_L += C2F_L[dir].ghost();
        std::cout << "L Span: " << D.span() << std::endl;
    }
    std::vector<Stencil<double>> C2F_H;
    Point span_H = Point::Zeros();
    for (int dir = 0; dir < DIM; dir++)
    {
        auto D = Stencil<double>::CellToFaceH(dir, Side::Lo, 5);
        //auto D = Stencil<double>::CellToFace(dir, Side::Hi, 5) * (1.0*Shift::Basis(dir, -1));
        C2F_H.push_back(D);
        span_H += C2F_H[dir].ghost();
        std::cout << "H Span: " << D.span() << std::endl;
    }
   
    int domainSize = 32;
    int numIter = 3;
    double err[numIter];
    double err_L[numIter];
    double err_H[numIter];
    for (int nn = 0; nn < numIter; nn++)
    {
        double dx = 1.0 / domainSize;
        Box domain = Box::Cube(domainSize);

        BoxData<double> phi_src(domain.grow(span));
        forallInPlace_p(f_wave_avg, phi_src, dx, k);
        h5.writePatch(dx, phi_src, "PHI_N%i", nn);

        BoxData<double, DIM> phi_soln(domain);
        forallInPlace_p(f_wave_avg_face, phi_soln, dx, k);
        h5.writePatch(dx, phi_soln, "PHI_SOLN_N%i", nn);
        
        BoxData<double, DIM> phi(domain);
        BoxData<double, DIM> phi_L(domain);
        BoxData<double, DIM> phi_H(domain);
        phi.setVal(7);
        phi_L.setVal(7);
        phi_H.setVal(7);

        for (int dir = 0; dir < DIM; dir++)
        {
            auto phi_i = slice(phi, dir);
            auto phi_L_i = slice(phi_L, dir);
            auto phi_H_i = slice(phi_H, dir);

            phi_i |= C2F[dir](phi_src);
            phi_L_i |= C2F_L[dir](phi_src);
            phi_H_i |= C2F_H[dir](phi_src);
        }

        BoxData<double, DIM> error(domain);
        BoxData<double, DIM> error_L(domain);
        BoxData<double, DIM> error_H(domain);
        phi.copyTo(error);
        phi_L.copyTo(error_L);
        phi_H.copyTo(error_H);
        error -= phi_soln;
        error_L -= phi_soln;
        error_H -= phi_soln;

        h5.writePatch(dx, phi,   "PHI_N%i", nn);
        h5.writePatch(dx, phi_L, "PHI_L_N%i", nn);
        h5.writePatch(dx, phi_H, "PHI_H_N%i", nn);
        h5.writePatch(dx, error,   "PHI_ERROR_N%i", nn);
        h5.writePatch(dx, error_L, "PHI_ERROR_L_N%i", nn);
        h5.writePatch(dx, error_H, "PHI_ERROR_H_N%i", nn);

        err[nn] = error.absMax();
        err_L[nn] = error_L.absMax();
        err_H[nn] = error_H.absMax();
        std::cout << "Error: " << err[nn] << std::endl;
        std::cout << "Error: " << err_L[nn] << std::endl;
        std::cout << "Error: " << err_H[nn] << std::endl;
        domainSize *= 2;
    }

    for (int ii = 1; ii < numIter; ii++)
    {
        std::cout << "Convergence Rate:     " << log(err[ii-1]/err[ii])/log(2.0) << std::endl;
        std::cout << "Convergence Rate (L): " << log(err_L[ii-1]/err_L[ii])/log(2.0) << std::endl;
        std::cout << "Convergence Rate (H): " << log(err_H[ii-1]/err_H[ii])/log(2.0) << std::endl;
    }
}
