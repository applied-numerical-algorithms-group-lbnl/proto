#include "Proto.H"
#include "InputParser.H"

using namespace Proto;

PROTO_KERNEL_START
void f_const_0(Point& a_pt, Var<double>& a_data, double a_dx, double a_value)
{
    a_data(0) = a_value;
}
PROTO_KERNEL_END(f_const_0, f_const);

PROTO_KERNEL_START
void f_sin_0(Point& a_pt, Var<double>& a_data, double a_dx)
{
    double x = a_pt[0]*a_dx + a_dx/2.0;
    double y = a_pt[1]*a_dx + a_dx/2.0;
    
    a_data(0) = sin(2.0*M_PI*(x + y));
}
PROTO_KERNEL_END(f_sin_0, f_sin);

PROTO_KERNEL_START
void f_tags_line_0(Point& a_pt, Var<short>& a_data, double a_dx, Point a_origin)
{
    Array<double, DIM> x;
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
PROTO_KERNEL_END(f_tags_line_0, f_tags_line);

PROTO_KERNEL_START
void f_tags_sphere_0(Point& a_pt, Var<short>& a_data, double a_dx, Point a_origin)
{
    Array<double, DIM> x;
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
PROTO_KERNEL_END(f_tags_sphere_0, f_tags_sphere);

PROTO_KERNEL_START
void f_tags_corner_0(Point& a_pt, Var<short>& a_data, double a_dx, Point a_corner)
{
    if (a_pt == a_corner)
    {
        a_data(0) = 1;
    } else {
        a_data(0) = 0;
    }
}
PROTO_KERNEL_END(f_tags_corner_0, f_tags_corner);

PROTO_KERNEL_START
void f_gaussian_0(Point& a_pt, Var<double>& a_data, double a_dx, Point a_origin, double a_sigma)
{
    Array<double, DIM> x;
    double rr = 0;
    for (int ii = 0; ii < DIM; ii++)
    {
        x[ii] = (a_pt[ii] - a_origin[ii])*a_dx + a_dx/2.0;
        rr += x[ii]*x[ii];
    }
   
    a_data(0) = exp(-rr / (2.0*a_sigma*a_sigma));
}
PROTO_KERNEL_END(f_gaussian_0, f_gaussian);

int main(int argc, char** argv)
{
#ifdef PR_MPI
    MPI_Init(&argc, &argv);
#endif
    using Proto::pout;
    HDF5Handler h5;
    
    int domainSize = 64;
    int boxSize = 32;
    double tagThreshold = 0.1;
    int refRatio = 2;
    int nestingDistance = 1;
    int TEST = -1;
    Array<bool, DIM> periodicity;
    for (int dir = 0; dir < DIM; dir++) { periodicity[dir] = true; }

    InputArgs args;
    
    args.add("TEST",        TEST);
    args.add("domainSize",      domainSize);
    args.add("boxSize",         boxSize);
    args.add("tagThreshold",    tagThreshold);
    args.add("refRatio",        refRatio);
    args.add("nestingDistance", nestingDistance);
    args.add("periodic_x",      periodicity[0]);
    args.add("periodic_y",      periodicity[1]);
    args.parse(argc, argv);
    args.print();

    if (procID() == 0)
    {
        if (TEST == -1)
        {
            std::cout << "Usage: ./meshRefine.exe -TEST <testNum>" << std::endl;
            std::cout << "Tests:" << std::endl;
            std::cout << "\tTest 0: Refine 2D Diagonal Line (Synthetic Tags)" << std::endl;
            std::cout << "\tTest 1: Refine ND Hollow Sphere (Synthetic Tags)" << std::endl;
            std::cout << "\tTest 2: Vorticity Distribution (Generate Tags From Data)" << std::endl;
            std::cout << "\tTest 10: Enforce Nesting In Bulk Domain" << std::endl;
            std::cout << "\tTest 11: Enforce Nesting At Periodic Boundaries" << std::endl;
        }
    }
#ifdef PR_MPI
    MPI_Bcast(&TEST, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    
    double L = 1.0;
    int bufferSize = 1;
    int one = 1;
    int zro = 0;
    Point refRatioV = Point::Ones(refRatio);
    for (int nn = 0; nn < 1; nn++)
    {
        double dx = L / domainSize;
        Point origin = Point::Ones(domainSize / 2);
       
        switch (TEST)
        {
            case 0:
            {
                Proto::pout() << "Running Test 0: 2D Diagonal Line" << std::endl;
                Box domain = Box::Cube(domainSize);

                Point tagBufferSize = Point::Ones(bufferSize);
                Point boxSizeVect = Point::Ones(boxSize);
                Point fineBoxSizeVect = Point::Ones(boxSize / 2); 

                Proto::pout() << "\tTag Buffer Size: " << tagBufferSize << std::endl;
                Proto::pout() << "\tBox Size (Coarse): " << boxSizeVect << std::endl;
                Proto::pout() << "\tBox Size (Fine): " << fineBoxSizeVect << std::endl;

                ProblemDomain problemDomain(domain, periodicity);
                DisjointBoxLayout layout(problemDomain, boxSizeVect);
                DisjointBoxLayout fineLayout(problemDomain.refine(refRatioV), fineBoxSizeVect);
                std::vector<DisjointBoxLayout> layouts;
                layouts.push_back(layout);
                layouts.push_back(fineLayout);
                std::vector<Point> refRatios;
                refRatios.push_back(refRatioV);

                AMRGrid grid(layouts, refRatios, 2);
                AMRData<double> testData(grid, Point::Zeros());
                testData[0].initialize(f_const, one, zro);
                testData[1].initialize(f_const, one, one);
                h5.writeAMRData({"data"}, 1.0, testData, "Grid_0");

                LevelTagData tags(layout, tagBufferSize);
                tags.initialize(f_tags_line, dx, origin);

                h5.writeLevel(tags, "Tags_0"); 
                AMRGrid::buffer(tags, tagBufferSize);
                h5.writeLevel(tags, "Tags_1"); 

                grid.regrid(tags, 0, fineBoxSizeVect);

                Point originPatch = origin / fineBoxSizeVect * refRatioV; 
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
                data[0].initialize(f_const, one, zro);
                data[1].initialize(f_const, one, one);

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
                DisjointBoxLayout fineLayout(problemDomain.refine(refRatioV), fineBoxSizeVect);
                std::vector<DisjointBoxLayout> layouts;
                layouts.push_back(layout);
                layouts.push_back(fineLayout);
                std::vector<Point> refRatios;
                refRatios.push_back(refRatioV);

                AMRGrid grid(layouts, refRatios, 2);
                AMRData<double> testData(grid, Point::Zeros());
                testData[0].initialize(f_const, one, zro);
                testData[1].initialize(f_const, one, one);
                h5.writeAMRData({"data"}, 1.0, testData, "Grid_0");

                LevelTagData tags(layout, tagBufferSize);
                tags.initialize(f_tags_sphere, dx, origin);

                h5.writeLevel(tags, "Tags_0"); 
                AMRGrid::buffer(tags, tagBufferSize);
                h5.writeLevel(tags, "Tags_1"); 

                grid.regrid(tags, 0, fineBoxSizeVect);

                Point originPatch = origin / fineBoxSizeVect * refRatioV; 
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
                data[0].initialize(f_const, one, zro);
                data[1].initialize(f_const, one, one);

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
                std::vector<Point> refRatios;
                refRatios.push_back(refRatioV);
                AMRGrid grid(layout, refRatios, 2);
                AMRData<double> data(grid, bufferSize);
                double sigma = 0.25;
                data.initialize(dx, f_gaussian, origin, sigma);
                h5.writeAMRData({"data"}, dx, data, "InputData");  
              
                // Compute Tags 
                LevelTagData tags;
                AMRGrid::computeTags(tags, data[0], tagBufferVect, tagThreshold);
                h5.writeLevel({"tags"}, dx, tags, "TagData");
               
                // Use Tags to Refine Mesh
                grid.regrid(tags, 0, fineBoxSizeVect);
                // this is just so we can see the output grid
                AMRData<double> gridData(grid, Point::Zeros());
                gridData[0].initialize(f_const, dx, zro);
                double dxFine = dx / refRatio;
                gridData[1].initialize(f_const, dxFine, one);
                h5.writeAMRData({"gridData"}, dx, gridData, "Grid");
                
                Point originPatch = origin / fineBoxSizeVect * refRatioV; 
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
                double dx_vect[3];
                dx_vect[0] = dx;
                dx_vect[1] = dx_vect[0] / refRatio;
                dx_vect[2] = dx_vect[1] / refRatio;
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
                std::vector<Point> refRatios;
                refRatios.push_back(refRatioV);
                refRatios.push_back(refRatioV);

                AMRGrid grid(layout, refRatios, 3);

                LevelTagData tags_0(layout, Point::Zeros());
                Point corner_0 = Point::Ones(boxSize) + Point::Basis(0, 1);
                tags_0.initialize(f_tags_corner, dx_vect[0], corner_0);
                h5.writeLevel(dx_vect[0], tags_0, "Tags_0");

                grid.regrid(tags_0, 0, boxSizeVect);

                LevelTagData tags_1(grid[1], Point::Zeros());
                Point corner_1 = (Point::Ones(boxSize) + Point::Basis(0, 1)) * refRatioV;
                tags_1.initialize(f_tags_corner, dx_vect[1], corner_1);
                h5.writeLevel(dx_vect[1], tags_1, "Tags_1");

                grid.regrid(tags_1, 1, fineBoxSizeVect);

                AMRData<double> data_before(grid, Point::Zeros());
                for (int ii = 0; ii < 3; ii++)
                {

                    Proto::pout() << "Initializing level " << ii << std::endl;
                    data_before[ii].initialize(f_const, dx_vect[ii], ii);
                    data_before[ii].layout().print();
                    h5.writeLevel(dx_vect[ii], data_before[ii], "Grid_L%i_0", ii);
                }
                pout() << "Writing AMR data " << std::endl;
                //h5.writeAMRData(dx_vect[0], data_before, "Grid_0");
                
                pout() << "Before enfore nesting" << std::endl;
                grid.enforceNesting(1, nestingDistance);
                pout() << "After enfore nesting" << std::endl;

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
                double dx_vect[3];
                dx_vect[0] = dx;
                dx_vect[1] = dx_vect[0] / refRatio;
                dx_vect[2] = dx_vect[1] / refRatio;

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
                std::vector<Point> refRatios;
                refRatios.push_back(refRatioV);
                refRatios.push_back(refRatioV);

                AMRGrid grid(layout, refRatios, 3);

                LevelTagData tags_0(layout, Point::Zeros());
                Point corner_0 = Point::Zeros();
                tags_0.initialize(f_tags_corner, dx_vect[0], corner_0);
                h5.writeLevel(dx_vect[0], tags_0, "Tags_0");

                grid.regrid(tags_0, 0, boxSizeVect);

                LevelTagData tags_1(grid[1], Point::Zeros());
                Point corner_1 = refRatioV;
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
        } // End switch (TEST)

    } // End for (n)
#ifdef PR_MPI
    MPI_Finalize();
#endif
}
