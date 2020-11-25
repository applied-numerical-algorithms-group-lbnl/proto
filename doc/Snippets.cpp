#include "./../include/Proto.H"

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
  }
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
  }
  PROTO_KERNEL_END(foo_p_temp, foo_p)
}

// globally defined functions are also valid:
PROTO_KERNEL_START
void bar_temp(Var<double>& arg_0, int arg_1)
{
  arg_0(0) = arg_1;
}
PROTO_KERNEL_END(bar_temp, bar)

// globally defined functions are also valid:
PROTO_KERNEL_START
void bar_p_temp(Point& a_p, Var<double>& arg_0, int arg_1)
{
  arg_0(0) = a_p[0]*arg_1;
}
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
    //  BoxIterator 
    //================
    //! [proto_boxiter]   
    Box iterBox(Point(1,2,3,4,5,6), Point(2,4,6,8,10,12));

    for (auto iter = iterBox.begin(); !iter.done(); ++iter)
    {
      // Do something for each Point in iterBox
      //    *iter is a Point instance
      //    iteration proceeds along lower axes fastest, starting with Bx::low()
      //    e.g. dimension 0 is fastest, 1 second fastest, etc. 
    }

    //if desired, the user can also iterate backwards:
    for (auto iter = iterBox.rbegin(); !iter.done(); --iter)
    {
      // Do things BACKWARDS
    }
    //! [proto_boxiter]   
    cout << "There is no test code for the BxIterator sippet." << endl;    
  } 
  {
    //====================================================================
    // CopyTo Example
    //================

    //! [proto_copyTo]
    // Defining inputs
    Box srcBox = Box::Cube(4);                        //[(0,..,0), (3,...,3)]
    Box destBox = Box::Cube(4).shift(Point::Ones());  //[(1,...,1), (4,...,4)]
    
    BoxData<double, 3, 3> Src(srcBox);  
    Src.setVal(7);                         //Initialize data as 7
    BoxData<double, 2, 2> Dest(destBox);   //Destination data is uninitialized

    Point copyShift = Point::Ones(2); //(2,...,2)
    Box srcCopyBox = Box::Cube(3);         //[(0,...,0), (2,...,2)]
    
    // Call copyTo
    // This statement copies the data of components {{1,2},{1,2},{0,0}} of Src within [(0,...,0), (2,...,2)]
    // into the components {{0,1},{0,1},{0,0}} of Dest within the Box [(2,...,2), (4,...,4)].
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
          if (a_dest(ii,jj) != 7 && srcCopyBox.shift(copyShift).contains(a_pt))
          {
            printf("CopyTo snippet is broken\n");
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
    Box srcBox = Box::Cube(4);                        //[(0,..,0), (3,...,3)]
    Box destBox = Box::Cube(4).shift(Point::Ones());  //[(1,...,1), (4,...,4)]
    
    BoxData<double, 3, 3> Src(srcBox);
    Src.setVal(7);
    BoxData<double, 2, 2> Dest(destBox);   //Destination data is uninitialized
    
    Point copyShift = Point::Ones(2); //(2,...,2)
    Box srcCopyBox = Box::Cube(3);         //[(0,...,0), (2,...,2)]
    
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
          if (a_dest(ii,jj) != 7 && srcCopyBox.shift(copyShift).contains(a_pt)) 
         {
           printf("LinearInOut snippet is broken\n");
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
    Box srcBox = Box::Cube(4);
    BoxData<double, 1, 2, 3> Src(srcBox);
    Src.setVal(17);
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
    Box srcBox = Box::Cube(4);
    BoxData<double, 1, 2, 3> Src(srcBox);
    Src.setVal(17);
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
    Box srcBox = Box::Cube(4); 
    BoxData<double,DIM+2> U(srcBox);
    U.setVal(1);
    const double gamma = 1.4;
      // Call forall to create and initialize a new BoxData<double, DIM+2> W 
    auto W1 = forall<double, DIM+2>(consToPrim, U, gamma);
      // Equivalent computation but without "auto"
    BoxData<double, DIM+2> W2 = forall<double, DIM+2>(consToPrim, U, gamma);
      // The Box on which the output is defined is the intersection of all the input BoxData. 
      //    In this case, there is only one BoxData - U - So the output will have the same Box as U.
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
    Box srcBox = Box::Cube(4);
    Box destBox = Box::Cube(3); // Defining the output domain
    BoxData<double,DIM+2> U(srcBox);
    U.setVal(1);
    const double gamma = 1.4;
      // Call forall to create and initialize a new BoxData<double, DIM+2> W 
    auto W1 = forall<double, DIM+2>(consToPrim, destBox, U, gamma);
      // Equivalent computation but without "auto"
    BoxData<double, DIM+2> W2 = forall<double, DIM+2>(consToPrim, destBox, U, gamma);
      // Here, the output will be defined on the input Box destBox.
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
    Box srcBox = Box::Cube(4);
    BoxData<double,DIM> X(srcBox);
    X.setVal(1);
    const double phase = M_PI/4.0;
      // Call forall to create and initialize a new BoxData<double, DIM+2> W 
    auto Y1 = forall_p<double>(sineFunc, X, phase);
      // Equivalent computation but without "auto"
    BoxData<double> Y2 = forall_p<double>(sineFunc, X, phase);
      // The Box on which the output is defined is the intersection of all the input BoxData. 
      //    In this case, there is only one BoxData - X - So the output will have the same Box as X.
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
    Box srcBox = Box::Cube(4);
    Box destBox = Box::Cube(3);
    BoxData<double,DIM> X(srcBox);
    X.setVal(1);
    const double phase = M_PI/4.0;
      // Call forall to create and initialize a new BoxData<double, DIM+2> W 
    auto Y1 = forall_p<double>(sineFunc, destBox, X, phase);
      // Equivalent computation but without "auto"
    BoxData<double> Y2 = forall_p<double>(sineFunc, destBox, X, phase);
      // Here, the output will be defined on the input Box destBox.
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
    Box srcBox = Box::Cube(4); 
    Box destBox = Box::Cube(3); 
    BoxData<double,DIM+2> U(srcBox);
    U.setVal(1);
    BoxData<double,DIM+2> W(destBox);
    const double gamma = 1.4;
      // Call forall to initialize W
    forallInPlace(consToPrim, W, U, gamma);
      // The Box on which the computation will be done is the intersection of all the input BoxData domains 
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
    Box srcBox = Box::Cube(4); 
    Box destBox = Box::Cube(3); 
    Box computeBox = Box::Cube(2); 
    BoxData<double,DIM+2> U(srcBox);
    U.setVal(1);
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
    Box srcBox = Box::Cube(4);
    Box destBox = Box::Cube(3);
    BoxData<double,DIM> X(srcBox);
    X.setVal(1);
    BoxData<double> Y(destBox);
    const double phase = M_PI/4.0;
      // Call forall to create and initialize Y 
    forallInPlace_p(sineFunc, Y, X, phase);
      // The Box on which the output is defined is the intersection of all the input BoxData. 
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
    Box srcBox = Box::Cube(4);
    Box destBox = Box::Cube(3);
    Box computeBox = Box::Cube(2);
    BoxData<double,DIM> X(srcBox);
    X.setVal(1);
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
  {
    //====================================================================
    // Stencil Build Example
    //=======================
    //! [proto_stencil_build]
    // This example builds and applies the standard 2*DIM + 1 Point Laplacian Stencil:
   
      // In 1D, this can be done in a single line: 
      //       stencil   =  1 x low side          - 2 x center         + 1 x high side 
    Stencil<double> L_1D = 1.0*Shift::Basis(0,-1) - 2.0*Shift::Zeros() + 1.0*Shift::Basis(0,1);
      // In arbitrary dimension:
    Stencil<double> L = (-2.0*DIM)*Shift::Zeros();
    for (int dir = 0; dir < DIM; dir++)
    {
      L += 1.0*Shift::Basis(dir,1) + 1.0*Shift::Basis(dir,-1);
    }
    //! [proto_stencil_build]
  }
  {
    //====================================================================
    // Stencil Apply Example
    //=======================
    //! [proto_stencil_apply]
      //  Build a Stencil from the built in Proto Stencil library:
    auto L = Stencil<double>::Laplacian(); //2nd order 2*DIM + 1 Point Laplacian
      //  Initialize some source data
      //  For this example, we define our input based on our desired ouput: an 8^DIM Point cube. 
    int rangeSize = 8;
    Box rangeBox = Box::Cube(rangeSize);  // [(0,...,0), (7,...,7)]
      //  The domainBox can be computed from the desired range, tanking ghost (halo) cells into account
    Box domainBox = L.domain(rangeBox);  // L needs 1 layer of ghost cells, so the domainBox is: [ (-1,...,-1), (8,....,8)]
    Box computeBox = rangeBox.grow(-1);
      //  Initialize the data using forall
    double dx = 0.5*M_PI/rangeSize; //assuming a range of PI/2
    BoxData<double> Src = forall_p<double>([=] PROTO_LAMBDA (Point& pt, Var<double>& src)
    {
      src(0) = sin(pt[0]*dx); //set Src = sin(x)
    }, domainBox);            //don't forget the domain Box!
    
    //  Apply L to Src
      //  Method 1: Build a new BoxData from the Stencil computation
    BoxData<double> Dest_0 = L(Src);
      //  When no Box input is specified, the output will be defined on the largest possible Box (e.g. rangeBox)
    BoxData<double> Dest_1 = L(Src, computeBox);
      //  When a VALID Box input is specified, the output will be restricted to that Box.
      //  *Usually* it makes more sense to let Stencil do its magic and range Box for you.
      
      //  Method 2: Apply Stencil computation output to an existing BoxData
    BoxData<double> Dest_3(rangeBox);
    BoxData<double> Dest_4(rangeBox);
    Dest_4.setVal(0);
      //  REPLACE the data in the destination with the output of L(Src)
    Dest_3 |= L(Src);
      //  ADD the output of L(Src) to the existing data in the destination
    Dest_4 += L(Src); 
      //  Both the ADD and REPLACE operations can specify a compute Box as well if desired
      //  Again, specifying the Box input is not recommended unless it is necessary

    //  WARNING: Note that 'auto' does NOT work for Stencil evaluation!!
    //    Compile errors involving 'LazyStencil' could be caused by inappropriate use of 'auto'.
    //! [proto_stencil_apply]
  }
  {
    //====================================================================
    // Stencil Source Refine Example
    //===============================
    //! [proto_stencil_average]
    //  This example Stencil computes a linear average from fine data source data onto a coarse grid
    Stencil<double> Avg;
    int refRatio = 2;
    Box offsetBox = Box::Cube(refRatio);
    for (auto iter = offsetBox.begin(); iter != offsetBox.end(); ++iter)
    {
      Avg += 1.0*Shift(*iter);
    }
    Avg *= (1.0/offsetBox.size());
    //  WARNING: Stencils with different src/dest refinement ratios cannot be added, subtracted, etc.
    //    When building a Stencil with non-trivial refinement, set the refinement ratio last.
    Avg.srcRatio() = Point::Ones(refRatio);

    int rangeSize = 8;
    Box indexBox = Box::Cube(rangeSize);          //[(0,...,0), (7,....,7)]
    Box domainBox = indexBox.refine(refRatio);   //[(0,...,0), (15,...,15)]
    auto Src = forall_p<double>([] PROTO_LAMBDA (Point& pt, Var<double>& src)
    {
      src(0) = 0.0;
      for (int ii = 0; ii < DIM; ii++)
      {
        src(0) += pt[ii];
      }
    },domainBox);

    //  The index box represents the points at which to compute the Stencil.
    //    When the source is refined by r, the data used to compute the solution at point p and p + 1
    //    In the index box is shifted by r.
    
    //    In this example, r = 2. In 1D, the average at Point (0) is computed by:
    //      avg(0) = (1.0*src(0) + 1.0*src(1))/2.
    //    At the next Point, the forumula is:
    //      avt(1) = (1.0*src(2) + 1.0*src(3))/2.
    //    Manifesting a shift in the source data by r = 2.

    BoxData<double> Dest_0 = Avg(Src,indexBox);
    //  OR
    BoxData<double> Dest_1 = Avg(Src); //Stencil automatically determines the largest possible Box for Dest_1, given the data available.
    //  The "|=" and "+=" operators can be used here as well with identical semantics.
    //  For an example illustrating the usage of the destination refinement ratio, see the documentation for InterpStencil

    //! [proto_stencil_average]
  }
  {
    //====================================================================
    // Stencil Dest Refine Example
    //=============================
    //! [proto_stencil_dest_refine]
    int refRatio = 2;
    //  This is a very simple stencil; it scales a Point by 7.
    Stencil<double> S = 7.0*Shift::Zeros();
    S.destRatio() = Point::Ones(2);
    S.destShift() = Point::Ones();
    
    int domainSize = 8;
    Box indexBox = Box::Cube(domainSize);           //[(0,...,0), (7,....,7)]
    Box rangeBox = indexBox.refine(refRatio);      //[(0,...,0), (15,...,15)]
    BoxData<double> Src(indexBox);
    Src.setVal(1);
    BoxData<double> Dest(rangeBox);
    Dest.setVal(0);

    Dest |= S(Src);
    //  In the case of a non-trivial source refinement ratio, the Stencil is shifted by srcRatio instead of by 1 for each
    //    consecutive Stencil evaluation. By contrast, when the destination refinement ratio is non-trivial, the Stencil only
    //    populates one in every destRefRatio Points, with the remainder of points being untouched by the operation. 

    //  In this particular example, at each Point in the indexBox, the following graphic represents the Stencil update (DIM = 2):
    //
    //    SOURCE                DESTINATION             DESTINATION
    //    +-------------+       +------+------+         +------+------+
    //    |             |       |      |      |         |      |      |
    //    |             |       | 0    | 0    |         | 0    | 7    |
    //    |      1      |   S   +------+------+   --->  +------+------+
    //    |             |       |      |      |         |      |      |
    //    |             |       | 0    | 0    |         | 0    | 0    |
    //    +-------------+       +------+------+         +------+------+
    //
    //! [proto_stencil_dest_refine]
  }
  {
    //====================================================================
    // InterpStencil Example
    //=======================
    //! [proto_stencil_interp]

    // Building the interpolation Stencil. This example uses piecewise constant interpolation
    int refRatio = 2;    
    InterpStencil<double> piecewiseConstant(refRatio);
    Box iterBox = Box::Cube(refRatio);  //[(0,...,0), (1,...,1)]
    for (auto destShift = iterBox.begin(); destShift != iterBox.end(); ++destShift)
    {
        // The InterpStencil is indexed into using the destShift value
        piecewiseConstant(*destShift) = 1.0*Shift::Zeros();
    }
    
    //This function sets all the shifts / ratios automatically and effectively makes the InterpStencil read-only
    //  Using InterpStencil::operator() on a closed InterpStencil will result in an error. Use InterpStencil::get(Point destShift) for read-only access
    piecewiseConstant.close();

    Box srcBox = Box::Cube(8);
    Box computeBox = srcBox;
    
    auto Src = forall_p<double>([] PROTO_LAMBDA (Point& pt, Var<double>& data)
    {
        data(0) = 0.0;
        for (int ii = 0; ii < DIM; ii++)
        {
            data(0) += pt[ii];
        }
    },srcBox);

    BoxData<double> Dest = piecewiseConstant(Src);
    //  In this particular example, at each Point in the indexBox, the following graphic represents the Stencil update (DIM = 2):
    //
    //    SOURCE                DESTINATION             DESTINATION
    //    +-------------+       +------+------+         +------+------+
    //    |             |       |      |      |         |      |      |
    //    |             |       | 0    | 0    |         | X    | X    |
    //    |      X      |   S   +------+------+   --->  +------+------+
    //    |             |       |      |      |         |      |      |
    //    |             |       | 0    | 0    |         | X    | X    |
    //    +-------------+       +------+------+         +------+------+
    //
    //  For all X in Src
    //! [proto_stencil_interp]
  }

  
cout << "If no message is saying otherwise, everything is working fine!" << endl; 
}
