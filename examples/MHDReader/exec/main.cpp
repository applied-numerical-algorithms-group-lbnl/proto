#include "Proto.H"
#include "MHDReader.H"

using namespace Proto;

int main(int argc, char** argv)
{
    #ifdef PR_MPI
    MPI_Init(&argc, &argv);
    #endif
    
    MHDReader reader;
    HDF5Handler h5;

    // Read data from "DATA.hdf5"
    BoxData<double, 8, HOST> data;
    std::vector<double> dtheta;
    reader.readData(data, "DATA");
    reader.readGeom(dtheta, "DATA");
    
    // Build grids
    Box domainBox(Point(180,360,180));
    Point boxSize = Point::Ones(90);
    std::array<bool, DIM> periodicity = {false, true, false}; //not read from HDF5 file
    ProblemDomain domain(domainBox, periodicity);
    DisjointBoxLayout layout(domain, boxSize);

    // Initialize data
    LevelBoxData<double, 8, HOST> state(layout, Point::Ones());
    state.setToZero(); // replace with initialization

    // Apply boundary condition
    for (auto index : layout)
    {
        auto& state_i = state[index];
        Point s(-1,0,0);            // shift in -x direction to put data in ghost region
        Box box_i = state_i.box().shift(s) & data.box(); //region of data to copy
        data.copyTo(state_i, box_i, s); 
    }

    h5.writeLevel(1, state, "STATE");
    //do something with dtheta...

    #ifdef PR_MPI
    MPI_Finalize();
    #endif
}
