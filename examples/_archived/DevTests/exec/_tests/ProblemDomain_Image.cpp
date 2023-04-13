#include "Proto.H"

using namespace Proto;

int main(int argc, char** argv)
{
    Box domainBox = Box::Cube(8);
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++)
    {
        if (dir == 0) { periodicity[dir] = false; }
        else { periodicity[dir] = true; }
    }
    ProblemDomain domain(domainBox, periodicity);
    
    Point p0 = Point::Ones();
    Point p1(0,-1,-1,-1,-1);
    Point p2(7,8,8,8,8,8,8);
    Point p3 = Point::Ones(-1);
    Point p4 = Point::Ones(8);

    std::cout << "Domain contains " << p0 << ": " << domain.contains(p0) << std::endl;
    std::cout << "Domain contains " << p1 << ": " << domain.contains(p1) << std::endl;
    std::cout << "Domain contains " << p2 << ": " << domain.contains(p2) << std::endl;
    std::cout << "Domain contains " << p3 << ": " << domain.contains(p3) << std::endl;
    std::cout << "Domain contains " << p4 << ": " << domain.contains(p4) << std::endl;

    Point q0 = domain.image(p0);
    Point q1 = domain.image(p1);
    Point q2 = domain.image(p2);

    std::cout << "Image of " << p0 << ": " << q0 << std::endl;
    std::cout << "Image of " << p1 << ": " << q1 << std::endl;
    std::cout << "Image of " << p2 << ": " << q2 << std::endl;
}
