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
    Point ghost = Point::Ones(4);
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    
    auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(domainSize);
    boxSizeVect[radialDir] = thickness;
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    auto map = CubedSphereShell::Map(layout, ghost);
    MBInterpOp iop;
    iop = CubedSphereShell::InterpOp<HOST>(layout, ghost);
    MBLevelBoxData<double, 1, HOST> UComp(layout, ghost);
    MBDisjointBoxLayout readLayout;
    h5.readMBLayout(readLayout, domain.graphPtr(), "UComp");
    layout.print();
    readLayout.print();
    h5.readMBLevel(UComp, "UComp");
    h5.writeMBLevel(map, UComp, "UCompPre");
    UComp.exchange();
    iop.apply(UComp,UComp);
    h5.writeMBLevel(map, UComp, "UCompPost");

#ifdef PR_MPI
    MPI_Finalize();
#endif
    return 0;
}
