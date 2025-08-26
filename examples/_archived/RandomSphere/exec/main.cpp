//#include <gtest/gtest.h>
#include "Proto.H"
//#include "RandomPointwiseFunctions.H"
//#include "MBLevelMap_Shear.H"
//#include "MBLevelMap_XPointRigid.H"
#include "MBLevelMap_CubeSphereShell.H"
#include "BoxOp_MBLaplace.H"

using namespace std;
#include <random>
#include "RandomPointwiseFunctions.H"

void GetCmdLineArgumenti(int argc, const char** argv, const char* name, int* rtn)
{
    size_t len = strlen(name);
    for(int i=1; i<argc; i+=2)
    {
        if(strcmp(argv[i]+1,name) ==0)
        {
            *rtn = atoi(argv[i+1]);
            std::cout<<name<<"="<<" "<<*rtn<<std::endl;
            break;
        }
    }
}
int main(int argc, char* argv[])
{
    HDF5Handler h5;
    int domainSize = 16;
    int boxSize = 16;
    int thickness = 1;
    bool cullRadialGhost = false;
    bool use2DFootprint = true;
    int numGhost = 5;
    int radialDir = CUBE_SPHERE_SHELL_RADIAL_COORD;
    Array<double, DIM> k{1,1,1,0,0,0};
    Array<double, DIM> offset{0,0,0,0,0,0};
    //offset += 0.1;
    Array<Point, DIM+1> dataGhost;
    Array<Point, DIM+1> errGhost;

    dataGhost.fill(Point::Ones(numGhost+2));
    dataGhost[0] = Point::Ones(numGhost);
    if (cullRadialGhost) { dataGhost[0][radialDir] = 0;}
    errGhost.fill(Point::Zeros());

    std::vector<Point> footprint;
    for (auto pi : Box::Kernel(3))
    {
        if (pi.abs().sum() <= 2)
        {
            if (use2DFootprint && (pi[radialDir] != 0)) { continue; }
            footprint.push_back(pi);
        }
    }
    int N = 1;
    double err[N];
    double errL1[N];
    for (int nn = 0; nn < N; nn++)
    {
        err[nn] = 0.0;
        errL1[nn] = 0.0;
        auto domain = buildCubeSphereShell(domainSize, thickness, radialDir);
        Point boxSizeVect = Point::Ones(boxSize);
        boxSizeVect[radialDir] = thickness;
        MBDisjointBoxLayout layout(domain, boxSizeVect);

        MBLevelMap_CubeSphereShell<HOST> map;
        map.define(layout, dataGhost);
        MBLevelMap_CubeSphereShellPolar<HOST> polarMaps[6];
        for (int bi = 0; bi < 6; bi++)
        {
            polarMaps[bi].define(layout, dataGhost, bi);
        }

        MBLevelBoxData<double, 1, HOST> hostSrc(layout, dataGhost);
        MBLevelBoxData<double, 1, HOST> hostDst(layout, dataGhost);
        MBLevelBoxData<double, 1, HOST> hostSln(layout, errGhost);
        MBLevelBoxData<double, 1, HOST> hostErr(layout, errGhost);

        hostSrc.setRandom(0,1);
        h5.writeMBLevel({"phi"}, map, hostSrc, "CubeSphereShell_Phi_Random_%i", nn);
    }
}
//#endif


