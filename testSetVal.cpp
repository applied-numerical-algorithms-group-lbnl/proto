#include "Proto.H"

using namespace Proto;

int main(int argc, char** argv)
{
    PR_TIMER_SETFILE("TestSetVal");
    int domainSize = 512;
    Box B = Box::Cube(domainSize);

    BoxData<double, DIM> data(B);
    {
        PR_TIME("oldSetVal");
        data.setVal(7.0);
    }
    {
        PR_TIME("newSetVal");
        std::fill(data.data(), data.data() + data.numValues(), 7.0);
    }
    PR_TIMER_REPORT();
    return 0;
}