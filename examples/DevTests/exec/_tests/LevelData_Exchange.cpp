#include "Proto.H"
#include "../../common/InputParser.H"
#include "func.H"

using namespace Proto;

void f_clearExterior(Point& a_pt, Var<double>& a_dst, Var<double>& a_src, Box a_box)
{
    if (a_box.contains(a_pt))
    {
        a_dst(0) = a_src(0);
    } else {
        a_dst(0) = 0;
    }
}



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

    LevelBoxData<double> srcData(layout, Point::Ones(ghostSize));
    LevelBoxData<double> dstData(layout, Point::Ones(ghostSize));
    LevelBoxData<double> errData(layout, Point::Ones(ghostSize));

    errData.setToZero(); 
    srcData.initialize(f_wave, dx, k); 
    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        auto& src_i = srcData[*iter];
        auto& dst_i = dstData[*iter];
        Box   box_i = iter.box();
        forallInPlace_p(f_clearExterior, dst_i, src_i, box_i);
    }

    h5.writeLevel(dx, srcData, "SRC");
    h5.writeLevel(dx, dstData, "DST_0");

    dstData.exchange();

    h5.writeLevel(dx, dstData, "DST_1");

    Box validDomainBox = domainBox;
    for (int dir = 0; dir < DIM; dir++)
    {
        if (periodicity[dir]) 
        {
            validDomainBox = validDomainBox.grow(Point::Basis(dir, ghostSize));
        }
    }
    Reduction<double, Abs> rxn;
    for (auto iter = layout.begin(); iter.ok(); ++iter)
    {
        auto& src_i = srcData[*iter];
        auto& dst_i = dstData[*iter];
        auto& err_i = errData[*iter];

        dst_i.copyTo(err_i);
        err_i -= src_i;
        forallInPlace_p(f_clearExterior, err_i, err_i, validDomainBox);
        
        err_i.absMax(rxn);
    }

    h5.writeLevel(dx, errData, "ERROR");

    if (rxn.fetch() < 1.0e-12) 
    {
        std::cout << "Error: " << rxn.fetch() << " | TEST PASSES." << std::endl;
    } else {
        std::cout << "Error: " << rxn.fetch() << " | TEST FAILS."  << std::endl;
    }
}
