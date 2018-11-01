#include "Proto.H"

using namespace Proto;
using namespace std;

//! [proto_forall_func]
// Valid funcion inputs to forall may be STATIC members of classes:
class Operator {

  // Pointwise function with no point dependence
  static void foo(Var<double, 3, 2>&  arg_0,
                  double              arg_1,  // plain-old-data can be passed by value
                  Var<bool>&          arg_2)  // any number of Var objects with different types / structures can be passed by reference
  {
    // if arg_2 == true at this Point...
    if (arg_2(0))
    {
      arg_0(1,1) = arg_1; // Access the (1,1,0) component at each point and set it to arg_1
    } else {
      arg_0(1,1) = -arg_1; // Access the (1,1,0) component at each point and se tit to -arg1
    }
  };

  // Pointwise function with point dependence
  static void foo_p(Point&              a_p,    // If the function depends on the point of evaluation, the Point must be the first argument
                    Var<double, 3, 2>&  arg_0,
                    Var<bool>&          arg_1)
  {
    if (arg_1(0))
    {
      for (int ii = 0; ii < DIM; ii++)
      {
        arg_0(1,1) += a_p[ii]; // Set the (1,1,0) component of arg_0 equal to the sum of the components of this Point
      }
    }
  };
};

// globally defined functions are also valid:
void bar(Var<double>& arg_0, int arg_1)
{
  arg_0(0) = arg_1;
};

// globally defined functions are also valid:
void bar_p(Point& a_p, Var<double>& arg_0, int arg_1)
{
  arg_0(0) = a_p[0]*arg_1;
};
//! [proto_forall_func

int main(int argc, char** argv)
{

  {
    //====================================================================
    // CopyTo Example
    //================

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
    forallInPlace_p([=] PROTO_LAMBDA (Point& a_pt, Var<double,2,2>& a_dest)
    {
      for (int ii = 0; ii < 2; ii++)
      {
        for (int jj = 0; jj < 2; jj++)
        {
          if ((a_dest(ii,jj) != 7) && srcCopyBox.shift(copyShift).contains(a_pt))
          {
            cout << "CopyTo snippet is broken" << endl;
            return;
          }
        }
      }
      //cout << "CopyTo snippet succeeded." << endl;
    }, Dest);
  }
  {
    //====================================================================
    // LinearIn / LinearOut Example
    //==============================
    
    //! [proto_linearInOut]
    Bx srcBox = Bx::Cube(4);                        //[(0,..,0), (3,...,3)]
    Bx destBox = Bx::Cube(4).shift(Point::Ones());  //[(1,...,1), (4,...,4)]
    
    BoxData<double, 3, 3> Src(srcBox, 7);  //Source data is initialized to 7
    BoxData<double, 2, 2> Dest(destBox);   //Destination data is uninitialized
    
    Point copyShift = Point::Ones(2); //(2,...,2)
    Bx srcCopyBox = Bx::Cube(3);         //[(0,...,0), (2,...,2)]
    
    double buffer[Src.box().size()*2*2];

    // Copy data from Src into the buffer
    Src.linearOut(buffer, srcCopyBox, {{1,2},{1,2},{}}); 
    
    // ... Operate on the buffer, send it in an MPI message, etc. ... 
    
    // Copy data from buffer into Dest
    Dest.linearIn(buffer, srcCopyBox.shift(copyShift), {{0,1},{0,1},{}}); 
    //! [proto_linearInOut]
    forallInPlace_p([=] PROTO_LAMBDA (Point& a_pt, Var<double,2,2>& a_dest)
    {
      for (int ii = 0; ii < 2; ii++)
      {
        for (int jj = 0; jj < 2; jj++)
        {
          if ((a_dest(ii,jj) != 7) && srcCopyBox.shift(copyShift).contains(a_pt))
          {
            cout << "LinearInOut snippet is broken" << endl;
            return;
          }
        }
      }
      //cout << "LinearInOut snippet succeeded." << endl;
    }, Dest);
  }
  {
    //====================================================================
    // Alias Example
    //================

    //! [proto_alias]
    Bx srcBox = Bx::Cube(4);
    BoxData<double, 1, 2, 3> Src(srcBox,17);
    // Alias is identical to Src and points to the same data. Changing alias will change Src.
    auto Alias = alias(Src);
    // shiftedAlias points to the same buffer as Src, but the domain is shifted by (1,...,1);
    //    (e.g. shiftedAlias[Point::Ones()] == Src[Point::Zeros] will return true.)
    auto shiftedAlias = alias(Src, Point::Ones());  //shiftedAlias points to the same data, but the associated domain
    //! [proto_alias]
    bool success = true;
    success &= (shiftedAlias.box() == srcBox.shift(Point::Ones()));
    for (auto iter = srcBox.begin(); iter != srcBox.end(); ++iter)
    {
      for (int ii = 0; ii < 1; ii++)
      for (int jj = 0; jj < 2; jj++)
      for (int kk = 0; kk < 3; kk++)
      {
        success &= (Alias.data(*iter,ii,jj,kk) == Src.data(*iter,ii,jj,kk));
        success &= (shiftedAlias.data((*iter + Point::Ones()),ii,jj,kk) == Src.data(*iter,ii,jj,kk));
      }
    }
    if (success) 
    {
      //cout << "Alias snippet was successful" << endl;
    }
    else 
    {
      cout << "Alias snippet is broken" << endl;
    }
  }
  {
    //====================================================================
    // Slice Example
    //================

    //! [proto_slice]
    Bx srcBox = Bx::Cube(4);
    BoxData<double, 1, 2, 3> Src(srcBox,17);
    //  Create an alias to the {1,1,0} component of Src
    //    Slice and Src are sharing data. slice[srcBox.low(),0,0,0] == Src[srcBox.low(),1,1,0] returns true;
    auto Slice = slice(Src, 0, 1);
    //! [proto_slice]

  }
  {
    //====================================================================
    // Forall Example
    //================
    

    Bx srcBox1 = Bx::Cube(4);                   //[(0,...,0), (3,...,3)]
    Bx srcBox2 = srcBox1.shift(Point::Ones());  //[(1,...,1), (4,...,4)]

    
  }
  cout << "If no message is saying otherwise, everything is working fine!" << endl; 
}
