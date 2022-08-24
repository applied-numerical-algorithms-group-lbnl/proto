#include "Proto.H"
#include "InputParser.H"
#include "PolarExchangeCopier.H"

using namespace Proto;

template<typename T, unsigned int C, MemType MEM>
    PROTO_KERNEL_START 
void f_initialize_(Point& a_pt, Var<T, C, MEM>& a_U, double a_dx)
{
    //Point x = a_pt + Point::Ones();
    double x = a_pt[0]*a_dx + a_dx/2.0;
    double y = a_pt[1]*a_dx + a_dx/2.0;
    a_U(0) = sin(2.0*M_PI*x)*sin(8.0*M_PI*y);
    //for (int cc = 0; cc < C; cc++)
    //{
    //    a_U(cc) = (cc+1) + 1e-2*x[0] + 1e-4*x[1];
    //}
#if DIM > 2
    double z = a_pt[2]*a_dx + a_dx/2.0;
    a_U(0) *= sin(4.0*M_PI*z);
    //for (int cc = 0; cc < C; cc++)
    //{
    //    a_U(cc) += 1e-6*x[2];
    //}
#endif
}
PROTO_KERNEL_END(f_initialize_, f_initialize)

int main(int argc, char** argv)
{

#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif

    // DEFAULT PARAMETERS
    int domainSize = 64;
    int boxSize = 16;
    
    // PARSE COMMAND LINE
    InputArgs args;
    args.add("domainSize",     domainSize);
    args.add("boxSize",        boxSize);
    args.parse(argc, argv);
    args.print();

    HDF5Handler h5;
    double dx = 1.0/domainSize;

    auto domain = Box::Cube(domainSize);
    array<bool,DIM> per;
    per.fill(true);
    ProblemDomain pd(domain,per);
    DisjointBoxLayout layout(pd,Point::Ones(boxSize));
    LevelBoxData<double> data(layout, Point::Ones(4));
    data.defineExchange<PolarExchangeCopier>(0,1);
    data.initialize(f_initialize, dx); 
    h5.writeLevel(dx, data, "DATA_0");
    data.exchange();
    h5.writeLevel(dx, data, "DATA_1");
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
