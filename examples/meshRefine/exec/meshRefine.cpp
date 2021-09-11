#include "Proto.H"

using namespace Proto;

void f_const (Point& a_pt, Var<double>& a_data, double a_dx, double a_value)
{
    a_data(0) = a_value;
}

void f_sin (Point& a_pt, Var<double>& a_data, double a_dx)
{
    double x = a_pt[0]*a_dx + a_dx/2.0;
    double y = a_pt[1]*a_dx + a_dx/2.0;
    
    a_data(0) = sin(2.0*M_PI*(x + y));
}

void f_tags_line (Point& a_pt, Var<char>& a_data, double a_dx, Point a_origin)
{
    std::array<double, DIM> x;
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = (a_pt[ii] - a_origin[ii])*a_dx + a_dx/2.0;
    }
    
    if (abs(x[0] - x[1]) < 2*a_dx) 
    {
        a_data(0) = 1;
    } else {
        a_data(0) = 0;
    }
}

void f_tags_sphere (Point& a_pt, Var<char>& a_data, double a_dx, Point a_origin)
{
    std::array<double, DIM> x;
    double r = 0;
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = (a_pt[ii] - a_origin[ii])*a_dx + a_dx/2.0;
        r += x[ii]*x[ii];
    }
    r = sqrt(r);
    double dr = 2*a_dx;
    
    if (r <= 0.3 + dr && r >= 0.3 - dr)
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
    int TEST_NUM = 0;
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
        if (argc > 3)
        {
            TEST_NUM = atoi(argv[3]);
        }
    }

#ifdef PR_MPI
    MPI_Bcast(&domainSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&boxSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&TEST_NUM, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    HDF5Handler h5;
    
    double L = 1.0;
    int bufferSize = 1;
    for (int nn = 0; nn < 1; nn++)
    {
        double dx = L / domainSize;
        Point origin = Point::Ones(domainSize / 2);
       
        switch (TEST_NUM)
        {
            case 0: {
            if (procID() == 0)
            {
                std::cout << "Running Test 0: 2D Diagonal Line" << std::endl;
            }
            Box domain = Box::Cube(domainSize);

            Point tagBufferSize = Point::Ones(bufferSize);
            Point boxSizeVect = Point::Ones(boxSize);
            Point fineBoxSizeVect = Point::Ones(boxSize / 2); 
            std::array<bool, DIM> periodicity;
            for (int ii = 0; ii < DIM; ii++)
            {
                if (ii == 0) { periodicity[ii] = true; }
                else {periodicity[ii] = false; }
            }
            
            if (procID() == 0)
            {
                std::cout << "\tTag Buffer Size: " << tagBufferSize << std::endl;
                std::cout << "\tBox Size (Coarse): " << boxSizeVect << std::endl;
                std::cout << "\tBox Size (Fine): " << fineBoxSizeVect << std::endl;
            }

            ProblemDomain problemDomain(domain, periodicity);
            DisjointBoxLayout layout(problemDomain, boxSizeVect);
            DisjointBoxLayout fineLayout(problemDomain.refine(Point::Ones(PR_AMR_REFRATIO)), fineBoxSizeVect);
            std::vector<DisjointBoxLayout> layouts;
            layouts.push_back(layout);
            layouts.push_back(fineLayout);

            AMRGrid grid(layouts, 2);
            AMRData<double> testData(grid, Point::Zeros());
            testData[0].initialize(f_const, 1, 0);
            testData[1].initialize(f_const, 1, 1);
            h5.writeAMRData({"data"}, 1.0, testData, "TestGrids");

            LevelTagData tags(layout, tagBufferSize);
            tags.initialize(f_tags_line, dx, origin);

            h5.writeLevel(tags, "Tags_0"); 
            AMRGrid::buffer(tags, tagBufferSize);
            h5.writeLevel(tags, "Tags_1"); 

            grid.regrid(tags, 0, fineBoxSizeVect);

            Point originPatch = origin / fineBoxSizeVect * PR_AMR_REFRATIO; 
            if (procID() == 0)
            {
                std::cout << "Created mesh is radially symmetric: " << grid[1].radialSymmetry(originPatch) << std::endl;
                for (int ii = 0; ii < DIM; ii++)
                {
                    std::cout << "Created mesh is symmetric in coordinate " << ii << ": " << grid[1].mirrorSymmetry(originPatch, ii) << std::endl;
                }
            }

            AMRData<double> data(grid, Point::Zeros());
            data.setToZero();

            h5.writeAMRData({"data"}, 1.0, data, "Grids");
            break;
            
            } default: {
            break;
            }
        } // End switch (TEST_NUM)

    } // End for (n)
    /*
       if (true)
       {
    // Create Finer Level Test
    Box domain = Box::Cube(domainSize);

    Point tagBufferSize = Point::Ones(bufferSize);
    ProblemDomain problemDomain(domain, true);
    DisjointBoxLayout layout(problemDomain, boxSizeVect);

    AMRGrid grid(layout, 2);
    LevelTagData tags(layout, tagBufferSize);
    tags.initialize(f_tags_line, dx, origin);

    h5.writeLevel(tags, "Tags_0"); 
    AMRGrid::buffer(tags, bufferSize);
    h5.writeLevel(tags, "Tags_1"); 

    grid.regrid(tags, 0);

    if (procID() == 0)
    {
    std::cout << "Created mesh is radially symmetric: " << grid[1].radialSymmetry(originPatch * PR_AMR_REFRATIO) << std::endl;
    for (int ii = 0; ii < DIM; ii++)
    {
    std::cout << "Created mesh is symmetric in coordinate " << ii << ": " << grid[1].mirrorSymmetry(originPatch * PR_AMR_REFRATIO, ii) << std::endl;
    }
    }

    AMRData<double> data(grid, Point::Zeros());
    data.setToZero();

    //LevelBoxData<double> data_0(grid[0], Point::Zeros());
    //LevelBoxData<double> data_1(grid[1], Point::Zeros());

    //data_0.setToZero();
    //data_1.setToZero();

    h5.writeLevel( 1.0, data[0], "Level_0");
    h5.writeLevel( 1.0/PR_AMR_REFRATIO, data[1], "Level_1");
    h5.writeAMRData({"data"}, 1.0, data, "AllLevels");

    } else {
    Point boxSizeVect = Point::Ones(boxSize);

    double dx_0 = dx;
    double dx_1 = dx_0 / PR_AMR_REFRATIO;
    double dx_2 = dx_1 / PR_AMR_REFRATIO;

    Box domain_0 = Box::Cube(domainSize);        
    ProblemDomain problemDomain_0(domain_0, true);
    DisjointBoxLayout layout_0(problemDomain_0, boxSizeVect);

    AMRGrid grid(layout_0, 3);

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
*/
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
