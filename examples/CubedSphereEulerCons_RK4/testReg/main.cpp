#include "Proto.H"
#include "Inputs_Parsing.H"
#include "BoxOp_EulerCubedSphere.H"
    template <typename T, MemType MEM>
BoxData<T,1,MEM> myRadialCoord(T a_dxi0, Box a_bx)
{
    int r_dir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    double r0 = CUBED_SPHERE_SHELL_R0;
    double r1 = CUBED_SPHERE_SHELL_R1;
    T dxi0 = a_dxi0;
    double dr = (r1-r0)*dxi0;
    BoxData<T,1,MEM> radius =
        forall_p<T,1,MEM>([]PROTO_LAMBDA(Point a_pt,
                    Var<double>& a_rad,
                    double a_r0,
                    double a_r1,
                    double a_dr,
                    double a_dxi0,
                    int a_rdir)
                {
                // T C_rad = 1.0;
                // T R_t = (a_r1 - a_r0) / (exp(C_rad) - 1.0);
                T etalow = a_pt[a_rdir]*a_dxi0 + .5*a_dxi0;
                // T rlow = a_r0 + R_t * (exp(C_rad * etalow) - 1.0);
                //a_rad(0) = rlow;
                a_rad(0) = a_r0 + etalow*(a_r1 - a_r0);
                },a_bx,r0,r1,dr,dxi0,r_dir);
    return radius;
}
int main(int argc, char *argv[])
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    ParseInputs::getInstance().parsenow(argc, argv);

    int domainSize = ParseInputs::get_domainSize();
    int thickness = ParseInputs::get_thickness();
    int boxSize_nonrad = ParseInputs::get_boxSize_nonrad();
    int boxSize_rad = ParseInputs::get_boxSize_rad();
    int max_iter = ParseInputs::get_max_iter();
    int write_cadence = ParseInputs::get_write_cadence();
    int convTestType = ParseInputs::get_convTestType();
    int init_condition_type = ParseInputs::get_init_condition_type();

    PR_TIMER_SETFILE(to_string(domainSize) + "_DIM" + to_string(DIM)  + "_CubeSphereTest.time.table");
    PR_TIMERS("MMBEuler");
    HDF5Handler h5;

    int lev = 0;
    typedef BoxOp_EulerCubedSphere<double, MBMap_CubedSphereShell, HOST> OP;
    int radialDir = CUBED_SPHERE_SHELL_RADIAL_COORD;
    auto domain = CubedSphereShell::Domain(domainSize, thickness, radialDir);
    Point boxSizeVect = Point::Ones(boxSize_nonrad);
    boxSizeVect[radialDir] = boxSize_rad;
    MBDisjointBoxLayout layout(domain, boxSizeVect);
    int countBox = 0;
    for (auto dit : layout)
    {
        countBox++;
    }
    std::cout << "proc_id: " << procID() << ";      num boxes: " << countBox << std::endl;

    auto map = CubedSphereShell::Map(layout, OP::ghost());
    Array<double, DIM> dx;
    dx[0] = 1.0/thickness;
    dx[1] = 1.0/domainSize;
    dx[2] = dx[1];

    MBBoundaryRegister<double,DIM,HOST> blockreg;
    blockreg.define(layout,1,Point::Zeros());
    blockreg.clear();

    // Fill registers. We compute DIM-component fluxes on each face, set to be the location
    // in Cartesian coordinates of the centers of the face.

    MBLevelBoxData<double, DIM, HOST> rhsTest(layout, Point::Zeros());
    MBLevelBoxData<double, DIM, HOST> levelFluxAv(layout, Point::Zeros());
    rhsTest.setVal(0.);
    levelFluxAv.setVal(0.);
    double half = 0.5;
    Stencil<double> av0 = 0.5*Shift(Point::Zeros()) + 0.5*Shift(Point::Basis(0));
    for (auto dit : layout)
    {
        Box bx = layout[dit];
        auto block = layout.block(dit);
        Box bxRadial = bx.extrude(Point::Ones() - Point::Basis(0));
        BoxData<double,1,HOST> radius = myRadialCoord<double,HOST>(dx[0],bxRadial);
        Array<BoxData<double,DIM,HOST>,DIM> fluxes;
        for (int dir = 1; dir < DIM ; dir++)
        {
            // Set fluxes to be the Cartesian spatial coordinate at cell faces.         
            Box bxFaces = bx.extrude(Point::Basis(dir));
            BoxData<double, DIM, HOST> XCart;           
            double offseta = half;
            double offsetb = half;            
            if (dir == 1) offseta = 0.;
            if (dir == 2) offsetb = 0.;
            XCart = forall_p<double,DIM,HOST>
                (f_cubedSphereMap3,bxFaces,radius,dx,offseta,offsetb,block);

            fluxes[dir].define(bxFaces);
            XCart.copyTo(fluxes[dir]);            
        }
        for (int dir = 1; dir < DIM; dir++)
        {
            Stencil<double> avdir = 0.5*Shift(Point::Zeros()) + 0.5*Shift(Point::Basis(dir));
            levelFluxAv[dit] += avdir(fluxes[dir],half);
        };

        blockreg.increment(fluxes,dit,1.0);
    }
    blockreg.exchange();
    blockreg.reflux(rhsTest);
    double maxReflux = rhsTest.absMax();
    double maxFlux = levelFluxAv.absMax();
    cout << "maxReflux/maxFlux = " << maxReflux / maxFlux << endl;
    std::string output_file = ParseInputs::get_data_file_prefix();
    h5.writeMBLevel({}, map, rhsTest, output_file);
    h5.writeMBLevel({}, map, levelFluxAv, "fluxAv");

    PR_TIMER_REPORT();
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
