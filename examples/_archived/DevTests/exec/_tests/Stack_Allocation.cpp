#include "Proto.H"

using namespace Proto;

int main(int argc, char** argv)
{
    Box b = Box::Cube(8);
    for (int ii = 0; ii < 3; ii++)
    {
        BoxData<double> data(b);
    }
}
