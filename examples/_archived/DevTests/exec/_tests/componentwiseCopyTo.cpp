#include "Proto.H"

using namespace Proto;

int main(int argc, char** argv)
{
    int boxSize = 4;
    int domainSize = 8;
    int k = 1;
    
    Box b = Box::Cube(4);
    BoxData<double, 3> A(b);
    BoxData<double, 2> B(b);
    BoxData<double, 2> C(b);
    

    A.setVal(1, b, 0);
    A.setVal(2, b, 1);
    A.setVal(3, b, 2);
    B.setVal(17);
    C.setVal(2, b, 0);
    C.setVal(3, b, 1);
    A.copyTo(B, b, {1,2}, Point::Zeros(), {0,1});

    A.printData();
    B.printData();

    C -= B;

    double error = C.absMax();
    if (error > 1e-12)
    {
        std::cout << "Error: " << error << " | TEST FAILED" << std::endl;
    } else {
        std::cout << "Error: " << error << " | TEST PASSED" << std::endl;
    }

    
}
