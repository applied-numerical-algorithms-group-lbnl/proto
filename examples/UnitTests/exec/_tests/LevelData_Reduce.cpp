#include "Proto.H"
#include "../../common/InputParser.H"
#include "func.H"

using namespace Proto;

int main(int argc, char** argv)
{
    HDF5Handler h5;

    int domainSize = 32;
    int boxSize = 8;
    int ghostSize = 1;
    std::array<bool, DIM> periodicity;

    InputArgs args;
    args.parse();
    args.set("domainSize", &domainSize);
    args.set("boxSize", &boxSize);
    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);

    double k = 1;
    double dx = 1.0 / domainSize;

    Point boxSizeV = Point::Ones(boxSize);
    Box domainBox = Box::Cube(domainSize);
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSizeV);

    LevelBoxData<double> data(layout, Point::Ones(ghostSize));
    data.initialize(f_wave, dx, k); 
    
    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        auto& data_i = data[*iter];
        Box   box_i = iter.box();
        forallInPlace_p(f_clearExterior, data_i, data_i, box_i);
    }

    h5.writeLevel(dx, data, "DATA_0");
    std::cout << "AbsMax: " << data.reduce<Abs>();
    std::cout << " | Max: " << data.reduce<Max>();
    std::cout << " | Min: " << data.reduce<Min>();
    std::cout << " | Sum: " << data.reduce<Sum>();
    std::cout << std::endl;

    data.exchange();
    h5.writeLevel(dx, data, "DATA_1");
    std::cout << "AbsMax: " << data.reduce<Abs>();
    std::cout << " | Max: " << data.reduce<Max>();
    std::cout << " | Min: " << data.reduce<Min>();
    std::cout << " | Sum: " << data.reduce<Sum>();
    std::cout << std::endl;
}
