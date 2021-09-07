#include "Proto.H"

using namespace Proto;

void f_sin (Point& a_pt, Var<double>& a_data, double a_dx)
{
    double x = a_pt[0]*a_dx + a_dx/2.0;
    double y = a_pt[1]*a_dx + a_dx/2.0;
    
    a_data(0) = sin(2.0*M_PI*(x + y));
}

void f_tags (Point& a_pt, Var<char>& a_data, double a_dx, Point a_origin)
{
    double x = (a_pt[0]-a_origin[0])*a_dx + a_dx/2.0;
    double y = (a_pt[1]-a_origin[0])*a_dx + a_dx/2.0;
    
    if (abs(x-y) < 3*a_dx)
    {
        a_data(0) = 1;
    } else {
        a_data(0) = 0;
    }
}

void f_tags_corner (Point& a_pt, Var<char>& a_data, double a_dx, Point a_corner)
{
    if (a_pt == a_corner)
    {
        a_data(0) = 1;
    } else {
        a_data(0) = 0;
    }
}

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    int domainSize = 8;
    int boxSize = 4;
    if (procID() == 0)
    {
        if (argc > 1)
        {
            domainSize = atoi(argv[1]);
        }
        if (argc > 2)
        {
            boxSize = atoi(argv[2]);
        }
    }

#ifdef PR_MPI
    MPI_Bcast(&domainSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&boxSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    HDF5Handler h5;

    double L = 1.0;
    double k = 2.0*M_PI;
    int refRatio = 2;

    for (int nn = 0; nn < 1; nn++)
    {
        double dx = L / domainSize;
        Point origin = Point::Ones(domainSize / 2);
        
        if (false)
        {
            // Create Finer Level Test
            Box domain = Box::Cube(domainSize);

            Point boxSizeVect = Point::Ones(boxSize);
            Point tagBufferSize = Point::Ones(2);
            ProblemDomain problemDomain(domain, true);
            DisjointBoxLayout layout(problemDomain, boxSizeVect);

            LevelTagData tags(layout, tagBufferSize);
            tags.initialize(f_tags, dx, origin);

            h5.writeLevel(tags, "Tags_0"); 
            AMRGrid::buffer(tags, 2);
            h5.writeLevel(tags, "Tags_1"); 

            AMRGrid grid(layout);
            grid.regrid(tags, 0);

            LevelBoxData<double> data_0(grid[0], Point::Zeros());
            LevelBoxData<double> data_1(grid[1], Point::Zeros());

            data_0.setToZero();
            data_1.setToZero();

            h5.writeLevel( 1.0, data_0, "Level_0");
            h5.writeLevel( 1.0/PR_AMR_REFRATIO, data_1, "Level_1");

        } else {
            Point boxSizeVect = Point::Ones(boxSize);
            
            double dx_0 = dx;
            double dx_1 = dx_0 / PR_AMR_REFRATIO;
            double dx_2 = dx_1 / PR_AMR_REFRATIO;
            
            Box domain_0 = Box::Cube(domainSize);        
            ProblemDomain problemDomain_0(domain_0, true);
            DisjointBoxLayout layout_0(problemDomain_0, boxSizeVect);
            
            AMRGrid grid(layout_0);
            
            LevelTagData tags_0(layout_0, Point::Zeros());
            Point corner_0 = Point::Ones(boxSize);
            tags_0.initialize(f_tags_corner, dx_0, corner_0);
            h5.writeLevel(dx_0, tags_0, "Tags_0");
            
            grid.regrid(tags_0, 0);
            
            LevelTagData tags_1(grid[1], Point::Zeros());
            Point corner_1 = Point::Ones(boxSize * PR_AMR_REFRATIO);
            tags_1.initialize(f_tags_corner, dx_1, corner_1);
            h5.writeLevel(dx_1, tags_1, "Tags_1");

            grid.regrid(tags_1, 1);
            
            LevelBoxData<double> data_L0_0(grid[0], Point::Zeros());
            LevelBoxData<double> data_L1_0(grid[1], Point::Zeros());
            LevelBoxData<double> data_L2_0(grid[2], Point::Zeros());

            data_L0_0.setToZero();
            data_L1_0.setToZero();
            data_L2_0.setToZero();

            h5.writeLevel(dx_0, data_L0_0, "Data_L0_0");
            h5.writeLevel(dx_1, data_L1_0, "Data_L1_0");
            h5.writeLevel(dx_2, data_L2_0, "Data_L2_0");

            grid.enforceNesting(1);

            LevelBoxData<double> data_L0_1(grid[0], Point::Zeros());
            LevelBoxData<double> data_L1_1(grid[1], Point::Zeros());
            LevelBoxData<double> data_L2_1(grid[2], Point::Zeros());

            data_L0_1.setToZero();
            data_L1_1.setToZero();
            data_L2_1.setToZero();

            h5.writeLevel(dx_0, data_L0_1, "Data_L0_1");
            h5.writeLevel(dx_1, data_L1_1, "Data_L1_1");
            h5.writeLevel(dx_2, data_L2_1, "Data_L2_1");



        }
        domainSize *= 2;
    }
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
