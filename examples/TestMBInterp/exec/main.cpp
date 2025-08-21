#include "TestMBInterp.H"
#include "Proto.H"

using namespace Proto;

int main(int argc, char* argv[])
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    HDF5Handler h5;
    int domainSize = 30;
    int thickness = 60;
    Point ghost = Point::Ones(7);
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    
    auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(domainSize);
    boxSizeVect[radialDir] = thickness;
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    auto map = CubedSphereShell::Map(layout, ghost);
    MBInterpOp iop = CubedSphereShell::InterpOp<HOST>(layout, ghost);
    Box B(Point(8,30,0), Point(8,36,36));
    for (auto pi : B)
    {
        MBDataPoint p(*layout.begin(), pi, layout);
        MBPointInterpOp& op = iop(p);
        std::cout << "Coefficients at Point " << pi << std::endl;
        for (auto ci : op.coefs())
        {
            std::cout << ci << ", ";
        }
        std::cout << std::endl;
        string fname = "INTERP_X" + std::to_string(pi[0]) + "_Y" + std::to_string(pi[1]) + "_Z" + std::to_string(pi[2]);
        op.writeFootprint(fname);
    }
    MBLevelBoxData<double, 1, HOST> UComp(layout, ghost);
    MBDisjointBoxLayout readLayout;
    h5.readMBLayout(readLayout, domain.graphPtr(), "UComp");
    layout.print();
    readLayout.print();
    h5.readMBLevel(UComp, "UComp");
    h5.writeMBLevel(map, UComp, "UCompPre");
    h5.writeMBLevelBoundsUnified({}, UComp, "UCompBounds_0");
    UComp.exchange();
    h5.writeMBLevelBoundsUnified({}, UComp, "UCompBounds_1");
    iop.apply(UComp,UComp);
    h5.writeMBLevel(map, UComp, "UCompPost");

#ifdef PR_MPI
    MPI_Finalize();
#endif
    return 0;
}
