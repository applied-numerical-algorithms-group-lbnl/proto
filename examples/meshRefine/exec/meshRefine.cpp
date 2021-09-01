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

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    int domainSize = 32;
    int boxSize = 8;
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


        domainSize *= 2;
    }
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
