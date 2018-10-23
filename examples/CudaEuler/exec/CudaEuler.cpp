#include "Proto.H"

using namespace Proto;

PROTO_KERNEL_START void foo_temp(Var<int>& data)
{
    data(0) = 1337;
}
PROTO_KERNEL_END(foo_temp,foo)

PROTO_KERNEL_START void bar_temp(Point& p, Var<int>& data)
{
    data(0) = p[0] + p[1];
}
PROTO_KERNEL_END(bar_temp,bar)

int main(int argc, char** argv)
{
    Bx B = Bx::Cube(4).shift(Point::Basis(0,2));
    BoxData<int> D0(B);
    cudaForall(foo,B,D0);
    D0.cudaSyncHost();
    D0.printData();
    
    BoxData<int> D1(B);
    cudaForall_p(bar,B,D1);
    D1.cudaSyncHost();
    D1.printData();
}
