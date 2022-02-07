#include "Proto.H"

using namespace Proto;

int main(int argc, char** argv)
{
    Box b = Box::Cube(4096);
    std::cout << "used: " << Stack<MEMTYPE_DEFAULT>::getStack().used() << std::endl;
    PR_STACK_ON;
    for (int ii = 0; ii < 10; ii++)
    {
        BoxData<double> data(b);
        std::cout << "allocating " << b.size()*sizeof(double) << std::endl;
        std::cout << "used: " << Stack<MEMTYPE_DEFAULT>::getStack().used() << std::endl;
    }
    PR_STACK_ON;
    for (int ii = 0; ii < 10; ii++)
    {
        BoxData<double> data(b);
        std::cout << "allocating " << b.size()*sizeof(double) << std::endl;
        std::cout << "used: " << Stack<MEMTYPE_DEFAULT>::getStack().used() << std::endl;
    }
    PR_STACK_OFF;
    std::cout << "used: " << Stack<MEMTYPE_DEFAULT>::getStack().used() << std::endl;
}
