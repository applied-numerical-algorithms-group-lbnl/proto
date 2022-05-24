#include "Proto.H"
#include "InputParser.H"

#define NUMCOMPS 1

using namespace Proto;

PROTO_KERNEL_START
void
f_fooF(const Point& a_pt, Var<int,NUMCOMPS, DEVICE>& a_data)
{
    Point x = a_pt + Point::Ones();
    for (int comp = 0; comp < NUMCOMPS; comp++)
    {
        a_data(comp) = (comp+1)*10 + x[0];
#if DIM > 1
        a_data(comp) = (comp+1)*100 + 10*x[0] + x[1];
#endif
#if DIM > 2
        a_data(comp) = (comp+1)*1000 + 100*x[0] + 10*x[1] + x[2];
#endif
    }
}
PROTO_KERNEL_END(f_fooF, f_foo);

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif
    
#ifdef PROTO_CUDA
    // SETUP
    HDF5Handler h5;

    int domainSize = 8;
    std::array<double, DIM> dx;
    dx.fill(1.0);

    Box domainBox = Box::Cube(domainSize);
    BoxData<int, NUMCOMPS, DEVICE> data_d(domainBox);

    forallInPlace_p(f_foo, data_d);

    h5.writePatch(dx, data_d, "DEVI_DATA");
#endif
    #ifdef PR_MPI
    MPI_Finalize();
    #endif
    return 0;
}

