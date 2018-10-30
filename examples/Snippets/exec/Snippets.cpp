#include "Proto.H"

using namespace Proto;

int main(int argc, char** argv)
{

  {
    //! [proto_copyTo]
    // Defining inputs
    Bx srcBox = Bx::Cube(4);                        //[(0,..,0), (3,...,3)]
    Bx destBox = Bx::Cube(4).shift(Point::Ones());  //[(1,...,1), (4,...,4)]
    
    BoxData<double, 3, 3> Src(srcBox, 7);  //Source data is initialized to 7
    BoxData<double, 2, 2> Dest(destBox);   //Destination data is uninitialized

    Point copyShift = Point::Ones(2); //(2,...,2)
    Bx srcCopyBox = Bx::Cube(3);         //[(0,...,0), (2,...,2)]
    
    // Call copyTo
    // This statement copies the data of components {{1,2},{1,2},{0,0}} of Src within [(0,...,0), (2,...,2)]
    // into the components {{0,1},{0,1},{0,0}} of Dest within the Bx [(2,...,2), (4,...,4)].
    Src.copyTo(Dest,              // Destination data
               srcCopyBox,        // Box to copy out of Src
               {{1,2},{1,2},{}},  // Components to copy out of Src defined by an in-place constructed CInterval
               copyShift,         // Amount to shift srcCopyBox by in the destination
               {{0,1},{0,1},{}}); // Components to copy into Dest defined by an in-place constructed CInterval
    //! [proto_copyTo]
  }
  
}
