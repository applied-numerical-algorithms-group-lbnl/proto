#include "Proto.H"
#include "InputParser.H"

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

void f_tags_line (Point& a_pt, Var<short>& a_data, double a_dx, Point a_origin)
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

void f_tags_sphere (Point& a_pt, Var<short>& a_data, double a_dx, Point a_origin)
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

void f_tags_corner (Point& a_pt, Var<short>& a_data, double a_dx, Point a_corner)
{
    if (a_pt == a_corner)
    {
        a_data(0) = 1;
    } else {
        a_data(0) = 0;
    }
}

void f_gaussian (Point& a_pt, Var<double>& a_data, double a_dx, Point a_origin, double a_sigma)
{
    std::array<double, DIM> x;
    double rr = 0;
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = (a_pt[ii] - a_origin[ii])*a_dx + a_dx/2.0;
        rr += x[ii]*x[ii];
    }
   
    a_data(0) = exp(-rr / (2.0*a_sigma*a_sigma));
}

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif

    HDF5Handler h5;
    InputArgs args;
    args.parse();
    args.print();
    
    int domainSize = 64;
    int boxSize = 32;
    int nestingDistance = 1;
    double tagThreshold = 0.1;
    std::array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    args.set("domainSize", &domainSize);
    args.set("boxSize", &boxSize);
    args.set("tagThreshold", &tagThreshold);
    args.set("nestingDistance", &nestingDistance);
    args.set("periodic_x", &periodicity[0]);
    args.set("periodic_y", &periodicity[1]);

    int TEST_NUM = -1;
    if (procID() == 0)
    {
        if (argc > 1)
        {
            TEST_NUM = atoi(argv[1]);
        }
        if (TEST_NUM == -1)
        {
            std::cout << "Usage: ./meshRefine.exe testNum" << std::endl;
            std::cout << "Tests:" << std::endl;
            std::cout << "\tTest 0: Refine 2D Diagonal Line (Synthetic Tags)" << std::endl;
            std::cout << "\tTest 1: Refine ND Hollow Sphere (Synthetic Tags)" << std::endl;
            std::cout << "\tTest 2: Vorticity Distribution (Generate Tags From Data)" << std::endl;
            std::cout << "\tTest 10: Enforce Nesting In Bulk Domain" << std::endl;
            std::cout << "\tTest 11: Enforce Nesting At Periodic Boundaries" << std::endl;
        }
    }

#ifdef PR_MPI
    MPI_Bcast(&TEST_NUM, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    
    double L = 1.0;
    int bufferSize = 1;
    for (int nn = 0; nn < 1; nn++)
    {
        double dx = L / domainSize;
        Point origin = Point::Ones(domainSize / 2);
       
        switch (TEST_NUM)
        {
            case 0:
            {
                pout() << "Running Test 0: 2D Diagonal Line" << std::endl;
                Box domain = Box::Cube(domainSize);

                Point tagBufferSize = Point::Ones(bufferSize);
                Point boxSizeVect = Point::Ones(boxSize);
                Point fineBoxSizeVect = Point::Ones(boxSize / 2); 

                pout() << "\tTag Buffer Size: " << tagBufferSize << std::endl;
                pout() << "\tBox Size (Coarse): " << boxSizeVect << std::endl;
                pout() << "\tBox Size (Fine): " << fineBoxSizeVect << std::endl;

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
                h5.writeAMRData({"data"}, 1.0, testData, "Grid_0");

                LevelTagData tags(layout, tagBufferSize);
                tags.initialize(f_tags_line, dx, origin);

                h5.writeLevel(tags, "Tags_0"); 
                AMRGrid::buffer(tags, tagBufferSize);
                h5.writeLevel(tags, "Tags_1"); 

                grid.regrid(tags, 0, fineBoxSizeVect);

                Point originPatch = origin / fineBoxSizeVect * PR_AMR_REFRATIO; 
                if (procID() == 0)
                {
                    std::cout << "Created mesh is radially symmetric: ";
                    std::cout << grid[1].radialSymmetry(originPatch) << std::endl;
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        std::cout << "Created mesh is symmetric in coordinate " << ii;
                        std::cout << ": " << grid[1].mirrorSymmetry(originPatch, ii) << std::endl;
                    }
                }

                AMRData<double> data(grid, Point::Zeros());
                data[0].initialize(f_const, 1, 0);
                data[1].initialize(f_const, 1, 1);

                h5.writeAMRData({"data"}, 1.0, data, "Grid_1");
                break;
            // ===========================================================
            // TEST 1
            } case 1: {
                if (procID() == 0)
                {
                    std::cout << "Running Test 1: ND Hollow Sphere" << std::endl;
                }
                Box domain = Box::Cube(domainSize);

                Point tagBufferSize = Point::Ones(bufferSize);
                Point boxSizeVect = Point::Ones(boxSize);
                Point fineBoxSizeVect = Point::Ones(boxSize); 

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
                h5.writeAMRData({"data"}, 1.0, testData, "Grid_0");

                LevelTagData tags(layout, tagBufferSize);
                tags.initialize(f_tags_sphere, dx, origin);

                h5.writeLevel(tags, "Tags_0"); 
                AMRGrid::buffer(tags, tagBufferSize);
                h5.writeLevel(tags, "Tags_1"); 

                grid.regrid(tags, 0, fineBoxSizeVect);

                Point originPatch = origin / fineBoxSizeVect * PR_AMR_REFRATIO; 
                if (procID() == 0)
                {
                    std::cout << "Created mesh is radially symmetric: ";
                    std::cout << grid[1].radialSymmetry(originPatch) << std::endl;
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        std::cout << "Created mesh is symmetric in coordinate " << ii;
                        std::cout << ": " << grid[1].mirrorSymmetry(originPatch, ii) << std::endl;
                    }
                }

                AMRData<double> data(grid, Point::Zeros());
                data[0].initialize(f_const, 1, 0);
                data[1].initialize(f_const, 1, 1);

                h5.writeAMRData({"data"}, 1.0, data, "Grid_1");
                break;
            // ===========================================================
            // TEST 2
            } case 2: {
                if (procID() == 0)
                {
                    std::cout << "Running Test 2: Vorticity Field" << std::endl;
                }
                Point bufferSize = Point::Ones();
                Box domain = Box::Cube(domainSize);
                Point boxSizeVect = Point::Ones(boxSize);
                Point fineBoxSizeVect = Point::Ones(boxSize / 2);
                Point tagBufferVect = Point::Zeros();
                if (procID() == 0)
                {
                    std::cout << "\tTag Buffer Size: " << tagBufferVect << std::endl;
                    std::cout << "\tBox Size (Coarse): " << boxSizeVect << std::endl;
                    std::cout << "\tBox Size (Fine): " << fineBoxSizeVect << std::endl;
                }

                // Initialize Data
                ProblemDomain problemDomain(domain, periodicity);
                DisjointBoxLayout layout(problemDomain, boxSizeVect);
                AMRGrid grid(layout, 2);
                AMRData<double> data(grid, bufferSize);
                data.initialize(dx, f_gaussian, origin, 0.25);
                h5.writeAMRData({"data"}, dx, data, "InputData");  
              
                // Compute Tags 
                LevelTagData tags;
                double tagThreshold = args.get("tagThreshold");
                AMRGrid::computeTags(tags, data[0], tagBufferVect, tagThreshold);
                h5.writeLevel({"tags"}, dx, tags, "TagData");
               
                // Use Tags to Refine Mesh
                grid.regrid(tags, 0, fineBoxSizeVect);
                // this is just so we can see the output grid
                AMRData<double> gridData(grid, Point::Zeros());
                gridData[0].initialize(f_const, dx, 0);
                gridData[1].initialize(f_const, dx / PR_AMR_REFRATIO, 1);
                h5.writeAMRData({"gridData"}, dx, gridData, "Grid");
                
                Point originPatch = origin / fineBoxSizeVect * PR_AMR_REFRATIO; 
                if (procID() == 0)
                {
                    std::cout << "Created mesh is radially symmetric: ";
                    std::cout << grid[1].radialSymmetry(originPatch) << std::endl;
                    for (int ii = 0; ii < DIM; ii++)
                    {
                        std::cout << "Created mesh is symmetric in coordinate " << ii;
                        std::cout << ": " << grid[1].mirrorSymmetry(originPatch, ii) << std::endl;
                    }
                }

                break;
            // ===========================================================
            // TEST 10
            } case 10: {
            
                if (procID() == 0)
                {
                    std::cout << "Running Test 10: Enforce Nesting In Bulk Domain" << std::endl;
                }
                int nestingDistance = args.get("nestingDistance");
                double dx_vect[3];
                dx_vect[0] = dx;
                dx_vect[1] = dx_vect[0] / PR_AMR_REFRATIO;
                dx_vect[2] = dx_vect[1] / PR_AMR_REFRATIO;
                Box domain = Box::Cube(domainSize);

                Point boxSizeVect = Point::Ones(boxSize);
                Point fineBoxSizeVect = Point::Ones(boxSize / 2); 

                if (procID() == 0)
                {
                    std::cout << "\tBox Size (Coarse): " << boxSizeVect << std::endl;
                    std::cout << "\tBox Size (Fine): " << fineBoxSizeVect << std::endl;
                }
                ProblemDomain problemDomain(domain, periodicity);
                DisjointBoxLayout layout(problemDomain, boxSizeVect);

                AMRGrid grid(layout, 3);

                LevelTagData tags_0(layout, Point::Zeros());
                Point corner_0 = Point::Ones(boxSize) + Point::Basis(0, 1);
                tags_0.initialize(f_tags_corner, dx_vect[0], corner_0);
                h5.writeLevel(dx_vect[0], tags_0, "Tags_0");

                grid.regrid(tags_0, 0, boxSizeVect);

                LevelTagData tags_1(grid[1], Point::Zeros());
                Point corner_1 = (Point::Ones(boxSize) + Point::Basis(0, 1)) * PR_AMR_REFRATIO;
                tags_1.initialize(f_tags_corner, dx_vect[1], corner_1);
                h5.writeLevel(dx_vect[1], tags_1, "Tags_1");

                grid.regrid(tags_1, 1, fineBoxSizeVect);

                AMRData<double> data_before(grid, Point::Zeros());
                for (int ii = 0; ii < 3; ii++)
                {
                    data_before[ii].initialize(f_const, dx_vect[ii], ii);
                }
                h5.writeAMRData(dx_vect[0], data_before, "Grid_0");

                grid.enforceNesting(1, nestingDistance);

                AMRData<double> data_after(grid, Point::Zeros());
                for (int ii = 0; ii < 3; ii++)
                {
                    data_after[ii].initialize(f_const, dx_vect[ii], ii);
                }
                h5.writeAMRData(dx_vect[0], data_after, "Grid_1");
                break;
            // ===========================================================
            // TEST 11
            } case 11: {
                if (procID() == 0)
                {
                    std::cout << "Running Test 10: Enforce Nesting At (Non)Periodic Boundaries" << std::endl;
                }
                int nestingDistance = args.get("nestingDistance");
                double dx_vect[3];
                dx_vect[0] = dx;
                dx_vect[1] = dx_vect[0] / PR_AMR_REFRATIO;
                dx_vect[2] = dx_vect[1] / PR_AMR_REFRATIO;

                Box domain = Box::Cube(domainSize);

                Point boxSizeVect = Point::Ones(boxSize);
                Point fineBoxSizeVect = Point::Ones(boxSize / 2); 

                if (procID() == 0)
                {
                    std::cout << "\tBox Size (Coarse): " << boxSizeVect << std::endl;
                    std::cout << "\tBox Size (Fine): " << fineBoxSizeVect << std::endl;
                }
                ProblemDomain problemDomain(domain, periodicity);
                DisjointBoxLayout layout(problemDomain, boxSizeVect);

                AMRGrid grid(layout, 3);

                LevelTagData tags_0(layout, Point::Zeros());
                //Point corner_0 = Point::Basis(1, boxSize / (2*PR_AMR_REFRATIO)) + Point::Ones();
                Point corner_0 = Point::Zeros();
                tags_0.initialize(f_tags_corner, dx_vect[0], corner_0);
                h5.writeLevel(dx_vect[0], tags_0, "Tags_0");

                grid.regrid(tags_0, 0, boxSizeVect);

                LevelTagData tags_1(grid[1], Point::Zeros());
                //Point corner_1 = (Point::Basis(1, boxSize / (2*PR_AMR_REFRATIO)) + Point::Ones()) * PR_AMR_REFRATIO;
                Point corner_1 = Point::Ones()*PR_AMR_REFRATIO;
                tags_1.initialize(f_tags_corner, dx_vect[1], corner_1);
                h5.writeLevel(dx_vect[1], tags_1, "Tags_1");

                grid.regrid(tags_1, 1, fineBoxSizeVect);

                AMRData<double> data_before(grid, Point::Zeros());
                for (int ii = 0; ii < 3; ii++)
                {
                    data_before[ii].initialize(f_const, dx_vect[ii], ii);
                }
                h5.writeAMRData(dx_vect[0], data_before, "Grid_0");

                grid.enforceNesting(1, nestingDistance);

                AMRData<double> data_after(grid, Point::Zeros());
                for (int ii = 0; ii < 3; ii++)
                {
                    data_after[ii].initialize(f_const, dx_vect[ii], ii);
                }
                h5.writeAMRData(dx_vect[0], data_after, "Grid_1");
                break;
            } default: {
                break;
            }
        } // End switch (TEST_NUM)

    } // End for (n)
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
