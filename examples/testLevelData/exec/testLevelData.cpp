#include "Proto.H"

using namespace Proto;

int main(int argc, char** argv)
{
    int domainSize = 8;
    Box domain = Box::Cube(domainSize);
    Point boxSize = Point::Ones(4);
    Point ghost = Point::Ones();
    std::array<bool, DIM> periodicity;
    for (int ii = 0; ii < DIM; ii++) { periodicity[ii] = false; }
    
    ProblemDomain problemDomain(domain, periodicity);
    DisjointBoxLayout layout(problemDomain, boxSize);

    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        std::cout << layout[*iter] << std::endl;
    }

    LevelBoxData<double> src(layout, ghost);
    LevelBoxData<double> dst(layout, ghost);

    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        src[*iter].setVal((*iter)+1, layout[*iter]);
        src[*iter].printData();
    }
    src.exchange();
    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        src[*iter].printData();
    }

}
