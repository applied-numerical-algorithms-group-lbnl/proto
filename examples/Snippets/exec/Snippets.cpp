#include "Proto.H"

using namespace Proto;
using namespace std;

//! [proto_forall_func]
// Valid funcion inputs to forall may be STATIC members of classes:
namespace Operator {

  // Pointwise function with no point dependence
  PROTO_KERNEL_START                          // necessary for use with GPU devices
  static void foo_temp(Var<double, 3, 2>&  arg_0,
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
  PROTO_KERNEL_END(foo_temp, foo)

  // Pointwise function with point dependence
  PROTO_KERNEL_START
  static void foo_p_temp(Point&              a_p,    // If the function depends on the point of evaluation, the Point must be the first argument
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
  PROTO_KERNEL_END(foo_p_temp, foo_p)
};

// globally defined functions are also valid:
PROTO_KERNEL_START
void bar_temp(Var<double>& arg_0, int arg_1)
{
  arg_0(0) = arg_1;
};
PROTO_KERNEL_END(bar_temp, bar)

// globally defined functions are also valid:
PROTO_KERNEL_START
void bar_p(Point& a_p, Var<double>& arg_0, int arg_1)
{
  arg_0(0) = a_p[0]*arg_1;
};
PROTO_KERNEL_END(bar_p_temp, bar_p)
//! [proto_forall_func]

//! [proto_forall_constoprim]
// This function is a physical example used in the Proto Euler example.
//    The input U contains the independent variables of the Euler equations in conservation form (e.g. rho*v_x)
//    The output W contains the corresponding primary variables (e.g. rho, v_x, v_y, etc.)
PROTO_KERNEL_START
void consToPrim_temp(Var<double,DIM+2>& W, 
                     const Var<double, DIM+2>& U,
                     double gamma)
{
  double rho = U(0);
  double v2 = 0.0;
  W(0) = rho;
  
  for (int i = 1; i <= DIM; i++)
    {
      double v;
      v = U(i) / rho;
      
      W(i) = v;
      v2 += v*v;
    }
  
  W(DIM+1) = (U(DIM+1) - .5 * rho * v2) * (gamma - 1.0);
}
PROTO_KERNEL_END(consToPrim_temp, consToPrim)
//! [proto_forall_constoprim]
//! [proto_forall_pointwise]
// Sample Function: Y = sine(x_0+phase) + sin(x_1 + phase) + ...
PROTO_KERNEL_START
void sineFunc_temp(Point& a_p,
                   Var<double>&     a_Y,
                   const Var<double,DIM>& a_X,
                   double           phase)
{
  a_Y(0) = 0.0;
  for (int ii = 0; ii < DIM; ii++)
  {
    a_Y(0) += sin(a_X(ii) + phase);
  }
}
PROTO_KERNEL_END(sineFunc_temp, sineFunc)
//! [proto_forall_pointwise]


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
    cout << "No test code written for slice" << endl;
  }
  {
    //====================================================================
    // Forall Example 1
    //==================
    //! [proto_forall_1]
      // Define inputs
    Bx srcBox = Bx::Cube(4); 
    BoxData<double,DIM+2> U(srcBox,1);
    const double gamma = 1.4;
      // Call forall to create and initialize a new BoxData<double, DIM+2> W 
    auto W1 = forall<double, DIM+2>(consToPrim, U, gamma);
      // Equivalent computation but without "auto"
    BoxData<double, DIM+2> W2 = forall<double, DIM+2>(consToPrim, U, gamma);
      // The Bx on which the output is defined is the intersection of all the input BoxData. 
      //    In this case, there is only one BoxData - U - So the output will have the same Bx as U.
      // The versions of forall that return a new BoxData have MANDATORY template arguments
      //    which must match those of the output BoxData.
    //! [proto_forall_1]
    cout << "No test code written for forall 1" << endl;
  }
  {
    //====================================================================
    // Forall Example 2
    //==================
    //! [proto_forall_2]
      // Define inputs
    Bx srcBox = Bx::Cube(4);
    Bx destBox = Bx::Cube(3); // Defining the output domain
    BoxData<double,DIM+2> U(srcBox,1);
    const double gamma = 1.4;
      // Call forall to create and initialize a new BoxData<double, DIM+2> W 
    auto W1 = forall<double, DIM+2>(consToPrim, destBox, U, gamma);
      // Equivalent computation but without "auto"
    BoxData<double, DIM+2> W2 = forall<double, DIM+2>(consToPrim, destBox, U, gamma);
      // Here, the output will be defined on the input Bx destBox.
      // The versions of forall that return a new BoxData have MANDATORY template arguments
      //    which must match those of the output BoxData.
    //! [proto_forall_2]
    
    cout << "No test code written for forall 2" << endl;
  }
  {
    //====================================================================
    // Forall Example 3
    //==================
    //! [proto_forall_3]
      // Define inputs
    Bx srcBox = Bx::Cube(4);
    BoxData<double,DIM> X(srcBox,1);
    const double phase = M_PI/4.0;
      // Call forall to create and initialize a new BoxData<double, DIM+2> W 
    auto Y1 = forall_p<double>(sineFunc, X, phase);
      // Equivalent computation but without "auto"
    BoxData<double> Y2 = forall_p<double>(sineFunc, X, phase);
      // The Bx on which the output is defined is the intersection of all the input BoxData. 
      //    In this case, there is only one BoxData - X - So the output will have the same Bx as X.
      // The versions of forall that return a new BoxData have MANDATORY template arguments
      //    which must match those of the output BoxData.
    //! [proto_forall_3]
    
    cout << "No test code written for forall 3" << endl;
  }
  {
    //====================================================================
    // Forall Example 4
    //==================
    //! [proto_forall_4]
      // Define inputs
    Bx srcBox = Bx::Cube(4);
    Bx destBox = Bx::Cube(3);
    BoxData<double,DIM> X(srcBox,1);
    const double phase = M_PI/4.0;
      // Call forall to create and initialize a new BoxData<double, DIM+2> W 
    auto Y1 = forall_p<double>(sineFunc, destBox, X, phase);
      // Equivalent computation but without "auto"
    BoxData<double> Y2 = forall_p<double>(sineFunc, destBox, X, phase);
      // Here, the output will be defined on the input Bx destBox.
      // The versions of forall that return a new BoxData have MANDATORY template arguments
      //    which must match those of the output BoxData.
    //! [proto_forall_4]
    
    cout << "No test code written for forall 4" << endl;
  }
  {
    //====================================================================
    // Forall Example 5
    //==================
    //! [proto_forall_5]
      // Define inputs
    Bx srcBox = Bx::Cube(4); 
    Bx destBox = Bx::Cube(3); 
    BoxData<double,DIM+2> U(srcBox,1);
    BoxData<double,DIM+2> W(destBox);
    const double gamma = 1.4;
      // Call forall to initialize W
    forallInPlace(consToPrim, W, U, gamma);
      // The Bx on which the computation will be done is the intersection of all the input BoxData domains 
      //    In this case, the domain will be srcBox & destBox = destBox
      // The "InPlace" versions of forall do not require template arguments
      //    In fact, they may not compile if the template arguments are provided, even if they are correct
    //! [proto_forall_5]
    cout << "No test code written for forall 5" << endl;
  }
  {
    //====================================================================
    // Forall Example 6
    //==================
    //! [proto_forall_6]
      // Define inputs
    Bx srcBox = Bx::Cube(4); 
    Bx destBox = Bx::Cube(3); 
    Bx computeBox = Bx::Cube(2); 
    BoxData<double,DIM+2> U(srcBox,1);
    BoxData<double,DIM+2> W(destBox);
    const double gamma = 1.4;
      // Call forall to initialize W
    forallInPlace(consToPrim, computeBox, W, U, gamma);
      // The computation is restricted to the Points in computeBox 
      // The "InPlace" versions of forall do not require template arguments
      //    In fact, they may not compile if the template arguments are provided, even if they are correct
    //! [proto_forall_6]
    cout << "No test code written for forall 6" << endl;
  }
  {
    //====================================================================
    // Forall Example 7
    //==================
    //! [proto_forall_7]
      // Define inputs
    Bx srcBox = Bx::Cube(4);
    Bx destBox = Bx::Cube(3);
    BoxData<double,DIM> X(srcBox,1);
    BoxData<double> Y(destBox);
    const double phase = M_PI/4.0;
      // Call forall to create and initialize Y 
    forallInPlace_p(sineFunc, Y, X, phase);
      // The Bx on which the output is defined is the intersection of all the input BoxData. 
      //    In this case, the computation is done on srcBox & destBox = destBox.
      // The "InPlace" versions of forall do not require template arguments
      //    In fact, they may not compile if the template arguments are provided, even if they are correct
    //! [proto_forall_7]
    
    cout << "No test code written for forall 7" << endl;
  }
  {
    //====================================================================
    // Forall Example 8
    //==================
    //! [proto_forall_8]
      // Define inputs
    Bx srcBox = Bx::Cube(4);
    Bx destBox = Bx::Cube(3);
    Bx computeBox = Bx::Cube(2);
    BoxData<double,DIM> X(srcBox,1);
    BoxData<double> Y(destBox);
    const double phase = M_PI/4.0;
      // Call forall to create and initialize Y 
    forallInPlace_p(sineFunc, computeBox, Y, X, phase);
      // The compution is restricted to computeBox
      // The "InPlace" versions of forall do not require template arguments
      //    In fact, they may not compile if the template arguments are provided, even if they are correct
    //! [proto_forall_8]
    
    cout << "No test code written for forall 8" << endl;
  }
  cout << "If no message is saying otherwise, everything is working fine!" << endl; 
}
