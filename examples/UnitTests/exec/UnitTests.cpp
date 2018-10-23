#include "Proto.H"
#include "UnitTest.H"
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

#define _USE_MATH_DEFINES

using namespace Proto;
using namespace std;
typedef Var<double,DIM> V;

//=================================================================================================
PROTO_KERNEL_START
void iotaFuncF(Point           & a_p,
               V               & a_X,
               double            a_h)
{
  for (int ii = 0; ii < DIM; ii++)
  {
    a_X(ii) = a_p[ii]*a_h + 0.5*a_h;
  }
}
PROTO_KERNEL_END(iotaFuncF,iotaFunc)

int main(int argc, char** argv)
{
  if (argc == 2)
  {
    VERBO = atoi(argv[1]);
  } else {
    VERBO = 0;
    cout << "You may call UnitTests.exe with a value of 1 or 2 to increase verbosity." << endl << endl;
  }

  cout << "What would you like to test?" << endl;
  cout << "\tTest Set 0: Run All Tests" << endl;
  cout << "\tTest Set 1: Point" << endl;
  cout << "\tTest Set 2: Box" << endl;
  cout << "\tTest Set 3: BoxData" << endl;
  cout << "\tTest Set 4: Stencil" << endl;
  cout << "\tTest Set 5: InterpStencil" << endl;
  int chosenTest;
  int numTests = 5;
  cin >> chosenTest;
  for (int ii = 1; ii <= numTests; ii++)
  {
    if (chosenTest == 0 || chosenTest == ii)
    {
      TEST = ii;
    }
    else
    {
      continue;
    }
    
    //***********************************
    //  POINT TESTS
    //***********************************
    BEGIN_TEST_SUITE(1,"Proto::Point");
    
    //===================================
    BEGIN_TEST("Default Construction");
    Point p;
    __OUT(2) cout << "Default constructed Point: " << p << endl; OUT__
    for (int ii = 0; ii < DIM; ii++)
    {
      UNIT_TEST((p[ii] == 0));
    }
    END_TEST();
    
    //===================================
    BEGIN_TEST("Array Construction");
    int v[DIM];
    for (int ii = 0; ii < DIM; ii++)
    {
        v[ii] = (ii+1)*17;
    }
    Point p(v);
    __OUT(2) {
      cout << "Building Array from [";
      for (int ii = 0; ii < DIM; ii ++)
      {
          cout << v[ii] << ", ";
      }
      cout << "]: " << p << endl;
    } OUT__
    for (int ii = 0; ii < DIM; ii++)
    {
      UNIT_TEST((p[ii] == v[ii]));
    }
    END_TEST();
    
    //===================================
    BEGIN_TEST("Variadic Construction");
    Point p(1,2,3,4,5,6);
    __OUT(2) {
      cout << "Building Point from [";
      for (int ii = 0; ii < DIM; ii ++)
      {
          cout << ii+1 << ", ";
      }
      cout << "]: " << p << endl;
    } OUT__
    for (int ii = 0; ii < DIM; ii++)
    {
      UNIT_TEST((p[ii] == ii+1));
    }
    END_TEST(); 
    
    //===================================
    BEGIN_TEST("Copy Construction");
    Point p(1,2,3,4,5,6);
    Point q(p);
    __OUT(2) {
        cout << "Constructing Point initial Point " << p <<": ";
        cout << q << endl;
    } OUT__
    for (int ii = 0; ii < DIM; ii++)
    {
        UNIT_TEST((p[ii] == q[ii]));
        UNIT_TEST(((&p[ii] != &q[ii])));
    }
    END_TEST();
    
    //===================================
    BEGIN_TEST("Static Methods");
    Point p0 = Point::Zeros();
    Point p1 = Point::Ones();
    Point p2 = Point::Ones(17);
    Point p3 = Point::Basis(0);
    Point p4 = Point::Basis(DIM-1,17);
    
    __OUT(2) {
        cout << "Point::Zeros(): " << p0 << endl; 
        cout << "Point::Ones(): " << p1 << endl; 
        cout << "Point::Ones(17): " << p2 << endl; 
        cout << "Point::Basis(0): " << p3 << endl; 
        cout << "Point::Basis("<<DIM-1<<", 17): " << p4 << endl; 
    } OUT__
    for (int ii = 0; ii < DIM; ii++)
    {
        UNIT_TEST((p0[ii] == 0));
        UNIT_TEST((p1[ii] == 1));
        UNIT_TEST((p2[ii] == 17));
        if (ii == 0){UNIT_TEST((p3[ii] == 1));}
        else {UNIT_TEST((p3[ii] == 0));}
        if (ii == DIM - 1){UNIT_TEST((p4[ii] == 17));}
        else {UNIT_TEST((p4[ii] == 0));}
    }
    END_TEST(); 

    //===================================
    BEGIN_TEST("Accessor Methods");
    Point p(1,2,3,4,5,6);
    for (int ii = 0; ii < DIM; ii++)
    {
        UNIT_TEST((p[ii] == ii+1));
    }
    const Point q(1,2,3,4,5,6);
    for (int ii = 0; ii < DIM; ii++)
    {
        UNIT_TEST((q[ii] == ii+1));
    }
    END_TEST();
   
    //===================================
    BEGIN_TEST("Algebraic Operators");
    Point p0(1,2,3,4,5,6);
    Point p1;
    int flip = -1;
    for (int ii = 0; ii < DIM; ii++)
    {
        p1[ii] = p0[ii]*17*flip;
        flip *= -1;
    }

    Point p2 = p0 + p1;
    __OUT(2) cout << p0 << " + " << p1 << " = " << p2 << endl; OUT__
    
    for (int ii = 0; ii < DIM; ii++)
    {
        UNIT_TEST((p2[ii] == (p0[ii] + p1[ii])));
    }
    
    p2 = p1 - p0;
    __OUT(2) cout << p1 << " - " << p0 << " = " << p2 << endl; OUT__
    
    for (int ii = 0; ii < DIM; ii++)
    {
        UNIT_TEST(((p2[ii] == (p1[ii] - p0[ii]))));
    }

    p2 = p1*p0;
    __OUT(2) cout << p1 << " * " << p0 << " = " << p2 << endl; OUT__

    for (int ii = 0; ii < DIM; ii++)
    {
        UNIT_TEST((p2[ii] == (p1[ii]*p0[ii])));
    }

    p2 = p1/p0;
    __OUT(2) cout << p1 << " / " << p0 << " = " << p2 << endl; OUT__

    for (int ii = 0; ii < DIM; ii++)
    {
        UNIT_TEST((p2[ii] == (p1[ii]/p0[ii])));
    }

    p2 = p1/17;
    __OUT(2) cout << p1 << " / 17 = " << p2 << endl; OUT__

    for (int ii = 0; ii < DIM; ii++)
    {
        UNIT_TEST((p2[ii] == (p1[ii]/17)));
    }

    p1 = Point::Ones(2); 
    p2 = p0 % p1;
    __OUT(2) cout << p0 << " % " << p1 << " = " << p2 << endl; OUT__

    for (int ii = 0; ii < DIM; ii++)
    {
        UNIT_TEST((p2[ii] == (p0[ii] % p1[ii])));
    }

    END_TEST();
     
    END_TEST_SUITE();
    
    //***********************************
    //  BX TESTS
    //***********************************
    BEGIN_TEST_SUITE(2,"Proto::Bx");
    //===================================
    BEGIN_TEST("Default Constructor");
    Bx B;
    __OUT(2) cout << "Default Constructed Bx: " << B << endl; OUT__
    UNIT_TEST((B.low() == Point::Zeros()));
    UNIT_TEST((B.high() == Point::Ones(-1)));
    UNIT_TEST((B.size() <= 0));
    END_TEST();

    //===================================
    BEGIN_TEST("Point Constructors");
    Point p1(1,2,3,4,5,6);
    Point p0 = (-1)*p1;
    Bx B0(p1,p0);
    Bx B1(p0,p1);
    Bx B2(p1);
    __OUT(2) {
        cout << "Bx(" << p1 << ", " << p0 << "): " << B0 << " ";
        cout << "size() = " << B0.size() << endl;
        cout << "Bx(" << p0 << ", " << p1 << "): " << B1 << " ";
        cout << "size() = " << B1.size() << endl;
        cout << "Bx(" << p1 << "): " << B2 << " ";
        cout << "size() = " << B2.size() << endl;
    } OUT__
    
    int s1 = 1; int s2 = 1;
    for (int ii = 0; ii < DIM; ii++)
    {
        s1 *= (p1[ii] - p0[ii] + 1);
        s2 *= (p1[ii] + 1);
    }
    UNIT_TEST((B0.size() <= 0));
    UNIT_TEST((B0.empty()));
    UNIT_TEST((B1.size() == s1));
    UNIT_TEST((!B1.empty()));
    UNIT_TEST((B2.size() == s2));
    UNIT_TEST((!B2.empty()));

    END_TEST();
    
    //===================================
    BEGIN_TEST("Copy Constructor");
    Bx B0(Point(1,2,3,4,5,6));
    Bx B1(B0);
    __OUT(2) cout << "Building Bx from " << B0 << ": " << B1 << endl; OUT__
    UNIT_TEST((B0 == B1));
    UNIT_TEST((&B0 != &B1)); 
    END_TEST();

    //===================================
    BEGIN_TEST("Static Methods");
    int size = 17;
    Bx B = Bx::Cube(size);
    if (VERBO > 1)
    {
        cout << "Bx::Cube(" << size << "): " << B << endl;
    }
    UNIT_TEST((B.size() == ipow<DIM>(size)));
    UNIT_TEST((B.low() == Point::Zeros()));
    UNIT_TEST((B.high() == Point::Ones(size-1)));
    END_TEST();

    //===================================
    BEGIN_TEST("Iteration and Indexing");
    Point max = Point(1,2,3,4,5,6);
    Point min = -1*max;
    Bx B(min,max);
    
    int index = 0;
    Point p = min;
    __OUT(2) cout << "Iterating through Box " << B << endl << endl; OUT__
    for (auto iter = B.begin(); iter != B.end(); ++iter, ++index)
    {
        __OUT(2) {cout << *iter << " ";} OUT__
        UNIT_TEST((*iter == p));
        UNIT_TEST((B.index(*iter) == index));
        ++p[0];
        for (int ii = 0; ii < DIM-1; ii++)
        {
            if (p[ii] > max[ii])
            {
                __OUT(2) {cout << endl;} OUT__
                p[ii] = min[ii];
                ++p[ii+1];
            }
        }
    }

    p = max;
    __OUT(2) cout << "Inverse Iterating through Box " << B << endl << endl; OUT__
    index = B.size()-1;
    for (auto iter = B.rbegin(); iter != B.rend(); --iter, --index)
    {
        __OUT(2) cout << *iter << " ";OUT__
        
        UNIT_TEST((*iter == p));
        UNIT_TEST((B.index(*iter) == index));
        --p[0];
        for (int ii = 0; ii < DIM-1; ii++)
        {
            if (p[ii] < min[ii])
            {
                __OUT(2) {cout << endl;} OUT__
                p[ii] = max[ii];
                --p[ii+1];
            }
        }
    }
    END_TEST();

    //===================================
    BEGIN_TEST("Shift");
    Bx b0 = Bx::Cube(16);
    Point s(1,-2,3,-4,5,-6);
    
    Bx b1 = b0.shift(s);
    if (VERBO > 1)
    {
        cout << "Shifting " << b0 << " by " << s << ": " << b1 << endl; 
    }
    UNIT_TEST(((b0.low() + s) == b1.low()));
    UNIT_TEST(((b0.high() + s) == b1.high()));
    UNIT_TEST((&b0 != &b1));
    UNIT_TEST((b0 == Bx::Cube(16))); //should be unchanged
    END_TEST();
    
    //===================================
    BEGIN_TEST("Grow");
    Bx b0 = Bx::Cube(16);
    Point s(1,-2,3,-4,5,-6);
    //Grow (Point)
    Bx b1 = b0.grow(s);
    if (VERBO > 1)
    {
        cout << "Growing " << b0 << " by " << s << ": " << b1 << endl;
    }
    UNIT_TEST(((b0.low() - s) == b1.low()));
    UNIT_TEST(((b0.high() + s) == b1.high()));
    UNIT_TEST((&b0 != &b1));
    UNIT_TEST((b0 == Bx::Cube(16))); //should be unchanged
    
    //Grow (scalar)
    b1 = b0.grow(3);
    if (VERBO > 1)
    {
        cout << "Growing " << b0 << " by " << 3 << ": " << b1 << endl;
    }
    UNIT_TEST(((b0.low() - 3) == b1.low()));
    UNIT_TEST(((b0.high() + 3) == b1.high()));
    UNIT_TEST((&b0 != &b1));
    UNIT_TEST((b0 == Bx::Cube(16))); //should be unchanged
    END_TEST();
    
    //===================================
    BEGIN_TEST("Coarsen"); 
    Point low = Point::Ones(-2);
    Point high = Point::Ones(3);
    Point r = Point::Ones() + Point::Basis(0); //(2,1,1,...,1)
    Bx b0 = Bx(low,high); 
    Bx b1 = b0.coarsen(2); 
    Bx b2 = b0.coarsen(3);
    Bx b3 = b0.coarsen(r);
    if (VERBO > 1)
    {
        cout << "Coarsen " << b0 << " by 2: " << b1 << endl;
        cout << "Coarsen " << b0 << " by 3: " << b2 << endl;
        cout << "Coarsen " << b0 << " by " << r << ": " << b3 << endl;
    }
    UNIT_TEST((b1.low() == b0.low()/2));
    UNIT_TEST((b1.high() == b0.high()/2));
    UNIT_TEST((b2.low() == b0.low()/3));
    UNIT_TEST((b2.high() == b0.high()/3));
    UNIT_TEST((b3.low()[0] == b0.low()[0]/2));
    UNIT_TEST((b3.high()[0] == b0.high()[0]/2));
    for (int ii = 1; ii < DIM; ii++)
    {
        UNIT_TEST((b3.low()[ii] == b0.low()[ii]));
        UNIT_TEST((b3.high()[ii] == b0.high()[ii]));
    }
    UNIT_TEST((&b0 != &b1));
    UNIT_TEST((&b0 != &b2));
    UNIT_TEST((&b0 != &b3));
    UNIT_TEST((b0 == Bx(low,high)));
    END_TEST();
    
    //===================================
    BEGIN_TEST("Refine"); 
    Point low = Point::Ones(-2);
    Point high = Point::Ones(3);
    Point r = Point::Ones() + Point::Basis(0); //(2,1,1,...,1)
    Bx b0 = Bx(low,high);
    Bx b1 = b0.refine(2);
    Bx b2 = b0.refine(3);
    Bx b3 = b0.refine(r);
    if (VERBO > 1)
    {
        cout << "Refine " << b0 << " by 2: " << b1 << endl;
        cout << "Refine " << b0 << " by 3: " << b2 << endl;
        cout << "Refine " << b0 << " by " << r << ": " << b3 << endl;
    }
    UNIT_TEST((b1.low() == b0.low()*2));
    UNIT_TEST((b1.high() == (b0.high()+Point::Ones())*2-Point::Ones()));
    UNIT_TEST((b2.low() == b0.low()*3));
    UNIT_TEST((b2.high() == (b0.high()+Point::Ones())*3-Point::Ones()));
    UNIT_TEST((b3.low()[0] == b0.low()[0]*2));
    UNIT_TEST((b3.high()[0] == (b0.high()[0]+1)*2 - 1));
    for (int ii = 1; ii < DIM; ii++)
    {
        UNIT_TEST((b3.low()[ii] == b0.low()[ii]));
        UNIT_TEST((b3.high()[ii] == b0.high()[ii]));
    }
    UNIT_TEST((&b0 != &b1));
    UNIT_TEST((&b0 != &b2));
    UNIT_TEST((&b0 != &b3));
    UNIT_TEST((b0 == Bx(low,high)));

    END_TEST();
    
    //===================================
    BEGIN_TEST("Flatten"); 
    Bx b0 = Bx::Cube(17);
    for (int ii = 0; ii < DIM; ii++)
    {
        Bx b1 = b0.flatten(ii);
        Bx b2 = b0.flatten(ii,true);
        if (VERBO > 1)
        {
            cout << "Down - flattening axis " << ii << " of " << b0;
            cout << ": " << b1 << endl;
            cout << "Up - flattening axis " << ii << " of " << b0;
            cout << ": " << b2 << endl;
        }
        UNIT_TEST((b1.low() == b0.low()));
        UNIT_TEST((b2.high() == b0.high()));
        for (int jj = 0; jj < DIM; jj++)
        {
            if (jj == ii)
            {
                UNIT_TEST((b1.high()[jj] == b1.low()[jj]));
                UNIT_TEST((b2.low()[jj] == b2.high()[jj]));
            } else {
                UNIT_TEST((b1.high()[jj] == b0.high()[jj]));
                UNIT_TEST((b2.low()[jj] == b0.low()[jj]));
            }
        }
    }
    END_TEST();

    //===================================
    BEGIN_TEST("Extrude"); 
    Bx b0 = Bx::Cube(17);
    for (int ii = 0; ii < DIM; ii++)
    {
        Bx b1 = b0.extrude(ii,3);
        Bx b2 = b0.extrude(ii,3,true);
        if (VERBO > 1)
        {
            cout << "Down - extruding axis " << ii << " of " << b0;
            cout << ": " << b1 << endl;
            cout << "Up - extruding axis " << ii << " of " << b0;
            cout << ": " << b2 << endl;
        }
        UNIT_TEST((b1.high() == b0.high()));
        UNIT_TEST((b2.low() == b0.low()));
        for (int jj = 0; jj < DIM; jj++)
        {
            if (jj == ii)
            {
                UNIT_TEST((b1.low()[jj] == b0.low()[jj]-3));
                UNIT_TEST((b2.high()[jj] == b0.high()[jj]+3));
            } else {
                UNIT_TEST((b1.low()[jj] == b0.low()[jj]));
                UNIT_TEST((b2.high()[jj] == b0.high()[jj]));
            }
        }
    }
    END_TEST();
    
    //===================================
    BEGIN_TEST("Mod"); 
    Point high(1,2,3,4,5,6);
    Point low = high*(-1);
    Bx B(low,high);
    
    Point p0 = Point::Zeros(); 
    Point p1(1,2,3,4,5,6);
    Point p2 = high+Point::Ones();
    Point p3 = low-Point::Ones();
    Point p4 = B.flatten(0).high() - Point::Basis(0);
    Point p5 = B.flatten(0,true).low() + Point::Basis(0);
    
    Point q0 = B.mod(p0);
    Point q1 = B.mod(p1);
    Point q2 = B.mod(p2);
    Point q3 = B.mod(p3);
    Point q4 = B.mod(p4);
    Point q5 = B.mod(p5);
    
    if (VERBO > 1)
    {
        cout << B << " mod " << p0 << " = " << q0 << endl;
        cout << B << " mod " << p1 << " = " << q1 << endl;
        cout << B << " mod " << p2 << " = " << q2 << endl;
        cout << B << " mod " << p3 << " = " << q3 << endl;
        cout << B << " mod " << p4 << " = " << q4 << endl;
        cout << B << " mod " << p5 << " = " << q5 << endl;
    }


    UNIT_TEST((q0 == p0));
    UNIT_TEST((q1 == p1));
    UNIT_TEST((q2 == low));
    UNIT_TEST((q3 == high));
    UNIT_TEST((q4 == high));
    UNIT_TEST((q5 == low));
    
    END_TEST();
    END_TEST_SUITE();
    
    //**********************************
    //  BOXDATA TESTS
    //**********************************
    BEGIN_TEST_SUITE(3,"Proto::BoxData");
    
    //===================================
    BEGIN_TEST("Default Constructor"); 
    BoxData<double,2,3> BD;
    __OUT(2) {
        cout << "Default Constructed BoxData<double,2,3>: " << endl;
        BD.print();
    } OUT__
    
    UNIT_TEST((BD.box() == Bx(Point::Ones(-1))));
    UNIT_TEST((BD.size() <= 0));
    END_TEST() 
   
    //===================================
    BEGIN_TEST("Box Constructor"); 
    Bx B = Bx(Point(1,2,3,4,5,6,7));
    BoxData<int, 3, 4, 5> BD(B);
    UNIT_TEST((BD.box() == B));
    int size = 1; 
    for (int ii = 0; ii < DIM; ii++)
    {
        size *= (ii+2);
    }
    UNIT_TEST((BD.size() == size*3*4*5));
    __OUT(2) {
        cout << "Uninitialized BoxData built from " << B << endl;
        BD.print();
    } OUT__
    END_TEST() 
   
    //===================================
    BEGIN_TEST("Initialization Constructor (And Accessor)"); 
    Bx B = Bx(Point(1,2,3,4,5,6,7));
    BoxData<int, 3, 4, 5> BD(B,1337);
    UNIT_TEST((BD.box() == B));
    int size = 1; 
    for (int ii = 0; ii < DIM; ii++)
    {
        size *= (ii+2);
    }
    UNIT_TEST((BD.size() == size*3*4*5));
    for (auto iter = B.begin(); iter != B.end(); ++iter)
    {
        for (int cc = 0; cc < 3; cc++)
        for (int dd = 0; dd < 4; dd++)
        for (int ee = 0; ee < 5; ee++)
        {
            UNIT_TEST((BD(*iter,cc,dd,ee) == 1337));
        }
    }
    __OUT(2) {
        cout << "BoxData built from " << B << " Initialized to 1337" << endl;
        BD.print();
    } OUT__
    END_TEST() 
   
    //===================================
    BEGIN_TEST("Move Constructor"); 

    #ifndef PROTO_MEM_CHECK
    OMIT_TEST("Move Constructor", "Compile with PROTO_MEM_CHECK=TRUE to run this test");
    #else
    FLUSH_CPY();
    Bx B = Bx(Point(1,2,3,4,5,6,7)*2);
    double dx = 0.1;
    auto X = forall_p<double,DIM>(iotaFunc,B,dx);
    UNIT_TEST((CPY.size() == 0));
    BoxData<double,DIM> Y(B,1337);
    FLUSH_CPY();
    Y = forall_p<double,DIM>(iotaFunc,B,dx);
    UNIT_TEST((CPY.size() == 0));
    for (auto iter = B.begin(); iter != B.end(); ++iter)
    {
        for (int ii = 0; ii < DIM; ii++)
        {
            UNIT_TEST((X(*iter,ii) == (*iter)[ii]*dx));
            UNIT_TEST((Y(*iter,ii) == ((*iter)[ii]*dx)));
        }
    }
    #endif
    END_TEST() 
   
    //===================================
    BEGIN_TEST("CInterval");
    __OUT(2) cout << "Testing General Constructor" << endl; OUT__
    CInterval I0(1,2,3,4,5,6);
    __OUT(2) {
      cout << "Building the CInterval I0(1,2,3,4,5,6): " << I0 << endl;
    } OUT__
    UNIT_TEST(I0.low(0) == 1);
    UNIT_TEST(I0.high(0) == 2);
    UNIT_TEST(I0.low(1) == 3);
    UNIT_TEST(I0.high(1) == 4);
    UNIT_TEST(I0.low(2) == 5);
    UNIT_TEST(I0.high(2) == 6);
    
    __OUT(2) cout << "Testing List Constructor" << endl; OUT__
    CInterval I1{{1,2},{3,4},{5,6}};
    CInterval I2{{},{3,4},{}};
    CInterval I3{1,2};
    __OUT(2) {
      cout << "Building the CInterval I1{{1,2},{3,4},{5,6}}: " << I1 << endl;
      cout << "Building the CInterval I2{{},{3,4},{}}: " << I2 << endl;
      cout << "Building the CInterval I3{1,2}: " << I3 << endl;
    } OUT__
    UNIT_TEST(I1.low(0) == 1);
    UNIT_TEST(I1.high(0) == 2);
    UNIT_TEST(I1.low(1) == 3);
    UNIT_TEST(I1.high(1) == 4);
    UNIT_TEST(I1.low(2) == 5);
    UNIT_TEST(I1.high(2) == 6);
    
    UNIT_TEST(I2.low(0) == 0);
    UNIT_TEST(I2.high(0) == 0);
    UNIT_TEST(I2.low(1) == 3);
    UNIT_TEST(I2.high(1) == 4);
    UNIT_TEST(I2.low(2) == 0);
    UNIT_TEST(I2.high(2) == 0);
   
    UNIT_TEST(I3.low(0) == 1);
    UNIT_TEST(I3.high(0) == 2);
    UNIT_TEST(I3.low(1) == 0);
    UNIT_TEST(I3.high(1) == 0);
    UNIT_TEST(I3.low(2) == 0);
    UNIT_TEST(I3.high(2) == 0);
  
    END_TEST();
   
    //===================================
    BEGIN_TEST("Copy"); 
    #ifndef PROTO_MEM_CHECK
    OMIT_TEST("Copy", "Compile with PROTO_MEM_CHECK=TRUE to run this test");
    #else
    Bx B = Bx(Point(1,2,3,4,5,6,7)*2);
    Point s = Point::Ones();
    Bx B1 = B.shift(s);
    Bx b = B.grow(-1);
    double dx = 0.1;
    auto X = forall_p<double,DIM>(iotaFunc, B, dx);
    
    BoxData<double,DIM> Y0(B,1337);
    BoxData<double,DIM> Y1(B1,1337);
    //CInterval<1> I0(0);
    //CInterval<1> I1(1);

    FLUSH_CPY();
    if (VERBO > 1)
    {
        cout << "Copying From: " << endl;
        X.printData();
        cout << "Copying To: " << endl;
        Y0.printData();
    }
    X.copyTo(Y0);
    if (VERBO > 1)
    {
        cout << "After Copy: " << endl;
        Y0.printData();
    }
    UNIT_TEST((CPY.size() == 1));
    UNIT_TEST((get<0>(CPY[0]) == &X));
    UNIT_TEST((get<1>(CPY[0]) == &Y0));
    for (auto iter = B.begin(); iter != B.end(); ++iter)
    {
        for (int cc = 0; cc < DIM; cc++)
        {
            UNIT_TEST((X(*iter,cc) == Y0(*iter,cc)));
            UNIT_TEST((&X(*iter,cc) != &Y0(*iter,cc)));
        }
    }
    
    FLUSH_CPY();
    if (VERBO > 1)
    {
        cout << "Copying From: " << endl;
        X.printData();
        cout << "Copying To: " << endl;
        Y1.printData();
    }
    X.copyTo(Y1);
    if (VERBO > 1)
    {
        cout << "After Copy: " << endl;
        Y1.printData();
    }
    UNIT_TEST((CPY.size() == 1));
    UNIT_TEST((get<0>(CPY[0]) == &X));
    UNIT_TEST((get<1>(CPY[0]) == &Y1));
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        if (X.box().contains(*iter))
        {
            for (int cc = 0; cc < DIM; cc++)
            {
                UNIT_TEST((X(*iter,cc) == Y1(*iter,cc)));
                UNIT_TEST((&X(*iter,cc) != &Y1(*iter,cc)));
            }
        } else {
            for (int cc = 0; cc < DIM; cc++)
            {
                UNIT_TEST((Y1(*iter,cc) == 1337));
            }
        }
    }
    Y1.setVal(1337);
    
    FLUSH_CPY();
    if (VERBO > 1)
    {

        cout << "Copying From: " << endl;
        X.printData();
        cout << "Copying To: " << endl;
        Y1.printData();
        cout << "In " << b << "Shifted by " << s << endl;
    }
    X.copyTo(Y1,b,s);
    if (VERBO > 1)
    {
        cout << "After Copy: " << endl;
        Y1.printData();
    }
    UNIT_TEST((CPY.size() == 1));
    UNIT_TEST((get<0>(CPY[0]) == &X));
    UNIT_TEST((get<1>(CPY[0]) == &Y1));
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        if (b.shift(s).contains(*iter))
        {
            for (int cc = 0; cc < DIM; cc++)
            {
                UNIT_TEST((X(*iter-s,cc) == Y1(*iter,cc)));
                UNIT_TEST((&X(*iter-s,cc) != &Y1(*iter,cc)));
            }
        } else {
            for (int cc = 0; cc < DIM; cc++)
            {
                UNIT_TEST((Y1(*iter,cc) == 1337));
            }
        }
    }
    Y1.setVal(1337);
    
    FLUSH_CPY();
    if (VERBO > 1)
    {
        cout << "Copying Component 0 From: " << endl;
        X.printData();
        cout << "Copying To Component 1: " << endl;
        Y1.printData();
        cout << "In " << b << "Shifted by " << s << endl;
    }
    X.copyTo(Y1,b,{0,0},s,{1,1});
    if (VERBO > 1)
    {
        cout << "After Copy: " << endl;
        Y1.printData();
    }
    UNIT_TEST((CPY.size() == 1));
    UNIT_TEST((get<0>(CPY[0]) == &X));
    UNIT_TEST((get<1>(CPY[0]) == &Y1));
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        if (b.shift(s).contains(*iter))
        {
            UNIT_TEST((X(*iter-s,0) == Y1(*iter,1)));
            UNIT_TEST((&X(*iter-s,0) != &Y1(*iter,1)));
            UNIT_TEST((Y1(*iter,0) == 1337));
        } else {
            for (int cc = 0; cc < DIM; cc++)
            {
                UNIT_TEST((Y1(*iter,cc) == 1337));
            }
        }
    }
    Y1.setVal(1337);
   
    //TODO: Need test with more than one tensor index
    
    #endif
    END_TEST() 
   
    //===================================
    BEGIN_TEST("Algebraic Operations"); 
    
    Bx B0 = Bx::Cube(4);
    Bx B1 = B0.shift(Point::Ones());
    double dx = 0.1;
    auto D0 = forall_p<double,DIM>(iotaFunc,B0,dx);
    BoxData<double,DIM> delta(B0,dx/2);
    D0 += delta;
    BoxData<double,DIM> D1(B1,17);

    if (VERBO > 1)
    {
        cout << "BoxData D0: " << endl;
        D0.printData();
        cout << "BoxData D1: " << endl;
        D1.printData();
    }

    D1 += D0;
    for (int cc = 0; cc < DIM; cc++)
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        if (B0.contains(*iter))
        {
            UNIT_TEST((D1(*iter,cc) == D0(*iter,cc) + 17));
        } else {
            UNIT_TEST((D1(*iter,cc) == 17));
        }
    }
    if (VERBO > 1)
    {
        cout << "D1 += D0" << endl;
        D1.printData();
    }
    D1.setVal(17);
    
    D1 -= D0;
    for (int cc = 0; cc < DIM; cc++)
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        if (B0.contains(*iter))
        {
            UNIT_TEST((D1(*iter,cc) == 17 - D0(*iter,cc)));
        } else {
            UNIT_TEST((D1(*iter,cc) == 17));
        }
    }
    if (VERBO > 1)
    {
        cout << "D1 -= D0" << endl;
        D1.printData();
    }
    D1.setVal(17);

    D1 *= D0;
    for (int cc = 0; cc < DIM; cc++)
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        if (B0.contains(*iter))
        {
            UNIT_TEST((D1(*iter,cc) == D0(*iter,cc) * 17));
        } else {
            UNIT_TEST((D1(*iter,cc) == 17));
        }
    }
    if (VERBO > 1)
    {
        cout << "D1 *= D0" << endl;
        D1.printData();
    }
    D1.setVal(17);
    
    D1 /= D0;
    for (int cc = 0; cc < DIM; cc++)
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        if (B0.contains(*iter))
        {
            UNIT_TEST((D1(*iter,cc) == 17 / D0(*iter,cc)));
        } else {
            UNIT_TEST((D1(*iter,cc) == 17));
        }
    }
    if (VERBO > 1)
    {
        cout << "D1 /= D0" << endl;
        D1.printData();
    }
    D1.setVal(17);

    END_TEST() 
   
    //===================================
    BEGIN_TEST("Reductions"); 
    
    Bx B = Bx::Cube(4).shift(Point::Basis(0,-2));
    double dx = 1;
    auto D = forall_p<double,DIM>(iotaFunc,B,dx);

    UNIT_TEST((D.max() == 3));
    UNIT_TEST((D.min() == -2));
    UNIT_TEST((D.absMax() == 3));
   
    for (int ii = 0; ii < DIM; ii++)
    {
        if (ii == 0)
        {
            UNIT_TEST((D.max(ii) == 1));
            UNIT_TEST((D.min(ii) == -2));
            UNIT_TEST((D.absMax(ii) == 2));
        } else {
            UNIT_TEST((D.max(ii) == 3));
            UNIT_TEST((D.min(ii) == 0));
            UNIT_TEST((D.absMax(ii) == 3));
        }
    }

    END_TEST() 
   
    //===================================
    BEGIN_TEST("Alias and Slice"); 
    #ifndef PROTO_MEM_CHECK
    OMIT_TEST("Alias and Slice","To run this test, compile with PROTO_MEM_CHECK=TRUE");
    #else
    Bx B0 = Bx::Cube(4).shift(Point::Basis(0,-2));
    Point shift = Point::Basis(0,-1);
    double dx = 0.1;
    auto D0 = forall_p<double,DIM>(iotaFunc,B0,dx);
    
    FLUSH_CPY();
    auto D1 = alias(D0,shift);
    UNIT_TEST((CPY.size() == 0));
    
    for (int ii = 0; ii < DIM; ii++)
    for (auto iter = B0.begin(); iter != B0.end(); ++iter)
    {
        UNIT_TEST((D1(*iter+shift,ii) == D0(*iter,ii)));
        UNIT_TEST((&D1(*iter+shift,ii) == &D0(*iter,ii)));
    }

    FLUSH_CPY();
    auto D2 = slice(D0,1);
    UNIT_TEST((CPY.size() == 0));
    
    for (auto iter = B0.begin(); iter != B0.end(); ++iter)
    {
        UNIT_TEST((D2(*iter,0) == D0(*iter,1)));
        UNIT_TEST((&D2(*iter,0) == &D0(*iter,1)));
    }

    #endif
    
    //TODO: Need slice tests with more than 1 tensor axis
    
    END_TEST() 
   
    //===================================
    BEGIN_TEST("Forall"); 
    #ifndef PROTO_MEM_CHECK
    OMIT_TEST("Forall","To run this test, please compile with PROTO_MEM_CHECK=TRUE");
    #else
    std::function<void(Var<double,DIM>&, const Var<double,DIM>&, const Var<int>&, Param<int>&)> foo = 
    [](Var<double,DIM>& a_v, const Var<double,DIM>& a_x, const Var<int>& a_c, Param<int>& a_n)
    {
        a_n()++;
        for (int ii = 0; ii < DIM; ii++)
        {
          a_v(ii) = a_x(ii) + a_c(0);
        }
    };
    
    const Bx B0 = Bx::Cube(5).shift(Point::Basis(0,-2));
    const Bx B1 = Bx::Cube(5);
    const Bx B2 = B0 & B1;
    const Bx b2 = Bx::Cube(2);
    double dx = 0.1;
    auto X = forall_p<double,DIM>(iotaFunc,B0,dx);
    BoxData<int> C(B1,17);
    int numPts = 0;

    if (VERBO > 1)
    {
        cout << "INPUTS:" << endl;
        cout << "X(x,y):" << endl;
        X.printData();
        cout << "C(c):" << endl;
        C.printData();
        cout << "Test problem:" << endl;
        cout << "OUT(x,y) = X(x,y) + C" << endl;
    }
    // forall
    //-------------------------------------------
    // with automatic Box
    if (VERBO > 0)
    {
      cout << endl << "Testing base forall" << endl << endl;
    }
    FLUSH_CPY();
    BoxData<double,DIM> D0 = forall<double,DIM>(foo,X,C,numPts);
    UNIT_TEST((CPY.size() == 0));
    
    UNIT_TEST((D0.box() == B2));
    UNIT_TEST((numPts == B2.size()));
    for (int dir = 0; dir < DIM; dir++)
    for (auto iter = B2.begin(); iter != B2.end(); ++iter)
    {
        UNIT_TEST((D0(*iter,dir) == C(*iter) + X(*iter,dir)));
    }
    if (VERBO > 1)
    {
        cout << "Output: " << b2 << endl;
        D0.printData();
    }
    
    // with supplied Box
    numPts = 0;
    if (VERBO > 0)
    {
      cout << endl << "Testing base forall with Box argument" << endl << endl;
    }
    FLUSH_CPY();
    BoxData<double,DIM> D1 = forall<double,DIM>(foo,b2,X,C,numPts);
    UNIT_TEST((CPY.size() == 0));
   
    UNIT_TEST((D1.box() == b2));
    UNIT_TEST((numPts == b2.size()));
    
    for (int dir = 0; dir < DIM; dir++)
    for (auto iter = b2.begin(); iter != b2.end(); ++iter)
    {
        UNIT_TEST((D1(*iter,dir) == C(*iter) + X(*iter,dir)));
    }

    if (VERBO > 1)
    {
        cout << "Output: " << b2 << endl;
        D1.printData();
    }

    //forallInPlace
    //-------------------------------------------
    // with automatic box
    BoxData<double,DIM> D2(B1,1337);
    BoxData<double,DIM> D3(B1,1337);

    if (VERBO > 0)
    {
      cout << endl << "Testing base in-place forall" << endl << endl;
    }
    numPts = 0;
    FLUSH_CPY();
    forallInPlace(foo,D2,X,C,numPts);
    UNIT_TEST((CPY.size() == 0));

    UNIT_TEST((numPts == B2.size()));
    for (int dir = 0; dir < DIM; dir++)
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        if (B2.contains(*iter))
        {
            UNIT_TEST((D2(*iter,dir) == X(*iter,dir) + C(*iter)));
        } else {
            UNIT_TEST((D2(*iter,dir) == 1337));
        }
    }
    if (VERBO > 1)
    {
      cout << "Output: " << endl;
      D2.printData();
    }
    // with box argument provided
    if (VERBO > 0)
    {
      cout << endl << "Testing in-place forall with box input" << endl << endl;
    }
    numPts = 0;
    FLUSH_CPY();
    forallInPlace(foo,b2,D3,X,C,numPts);
    UNIT_TEST((CPY.size() == 0));
    
    UNIT_TEST((numPts == b2.size()));

    for (int dir = 0; dir < DIM; dir++)
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        if (b2.contains(*iter))
        {
            UNIT_TEST((D3(*iter,dir) == X(*iter,dir) + C(*iter)));
        } else {
            UNIT_TEST((D3(*iter,dir) == 1337));
        }
    }
    if (VERBO > 1)
    {
      cout << "Output: " << endl;
      D3.printData();
    }
    #endif
    END_TEST();
    
    //===================================
    BEGIN_TEST("Forall_p");
    #ifndef PROTO_MEM_CHECK
    OMIT_TEST("Forall_p","To run this test, please compile with PROTO_MEM_CHECK=TRUE");
    #else
    
    std::function<void(Point&, Var<double>&, const Var<double>&)> square = 
    [](Point& a_p, Var<double>& a_v, const Var<double>& a_c)
    {
        double x = a_p[0];
        a_v(0) = x*x + a_c(0);
    };
    Bx B0 = Bx::Cube(8);
    Bx B1 = Bx::Cube(8).shift(Point::Basis(0,-1));
    Bx B2 = B0 & B1;
    Bx b2 = B2.grow(-1);
    
    BoxData<double> C(B0,0.17);
    
    if (VERBO > 1)
    {
        cout << "Test Problem: y = x^2 + 0.17" << endl;
    }
    if (VERBO > 0)
    {
      cout << endl << "Testing base forall_p" << endl << endl;
    }
    
    FLUSH_CPY();
    auto D0 = forall_p<double>(square,C);
    UNIT_TEST((CPY.size() == 0));
    
    UNIT_TEST((D0.box() == B0));
    
    for (auto iter = B0.begin(); iter != B0.end(); ++iter)
    {
        Point p = *iter;
        UNIT_TEST((D0(*iter) == p[0]*p[0]+0.17));
    }
    if (VERBO > 1)
    {
        cout << "Output: " << endl;
        D0.printData();
    }
    
    if (VERBO > 0)
    {
      cout << endl << "Testing base forall_p with box input" << endl << endl;
    }
    
    FLUSH_CPY();
    auto D3 = forall_p<double>(square,b2,C);
    UNIT_TEST((CPY.size() == 0));
    
    UNIT_TEST((D3.box() == b2));
    
    for (auto iter = b2.begin(); iter != b2.end(); ++iter)
    {
        Point p = *iter;
        UNIT_TEST((D3(*iter) == p[0]*p[0]+0.17));
    }
    if (VERBO > 1)
    {
        cout << "Output: " << endl;
        D3.printData();
    }

    BoxData<double> D1(B1,1337);
    BoxData<double> D2(B1,1337);
    if (VERBO > 0)
    {
      cout << endl << "Testing in place forall_p" << endl << endl;
    }
    FLUSH_CPY();
    forallInPlace_p(square,D1,C);
    UNIT_TEST((CPY.size() == 0));
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        Point p = *iter;
        if (B2.contains(p))
        {
            UNIT_TEST((D1(*iter) == p[0]*p[0] + 0.17));
        } else {
            UNIT_TEST((D1(*iter) == 1337));
        }
    }
    if (VERBO > 1)
    {
        cout << "Output: " << endl;
        D1.printData();
    }
    if (VERBO > 0)
    {
      cout << endl << "Testing in place forall_p with box input" << endl << endl;
    }
    
    FLUSH_CPY();
    forallInPlace_p(square,b2,D2,C);
    UNIT_TEST((CPY.size() == 0));
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        Point p = *iter;
        if (b2.contains(p))
        {
            UNIT_TEST((D2(*iter) == p[0]*p[0] + 0.17));
        } else {
            UNIT_TEST((D2(*iter) == 1337));
        }
    }
    if (VERBO > 1)
    {
        cout << "Output: " << endl;
        D2.printData();
    }
    #endif
    END_TEST();
     
    END_TEST_SUITE();
    
    //*********************************
    //  STENCIL TESTS
    //**********************************
    
    BEGIN_TEST_SUITE(4,"Proto::Stencil");
    
    //===================================
    BEGIN_TEST("Default Constructor"); 
    Stencil<double> S;
    Bx B = Bx::Cube(8);
    if (VERBO > 1)
    {
        cout << "Default Stencil: " << endl;
        S.print();
    }
    std::vector<Point> offsets = S.offsets(); 
    std::vector<double> coefs = S.coefs();
    UNIT_TEST((offsets.size() == 0));
    UNIT_TEST((coefs.size() == 0));
    UNIT_TEST((S.span() == Bx(Point::Zeros(), Point::Ones(-1)))); 
    END_TEST();
    
    //===================================
    BEGIN_TEST("Shift Constructor + Stencil Addition"); 
    Stencil<double> S = 1.0*Shift::Zeros() + 
                      + 5.0*Shift::Zeros() + 2.0*Shift::Ones()
                      + 3.0*Shift::Basis(0,2) + 4.0*Shift(6,5,4,3,2,1);
    if (VERBO > 1)
    {
        cout << "Building the Stencil: " << endl;
        cout << 1.0 << " * " << Point::Zeros() << " + ";
        cout << 5.0 << " * " << Point::Zeros() << " + ";
        cout << 2.0 << " * " << Point::Ones() << " + ";
        cout << 3.0 << " * " << Point::Basis(0,2) << " + ";
        cout << 4.0 << " * " << Point(6,5,4,3,2,1) << endl;
        S.print();
    }
    std::vector<Point> offsets = S.offsets(); 
    std::vector<double> coefs = S.coefs();
    
    UNIT_TEST((offsets.size() == 4));
    UNIT_TEST((coefs.size() == 4));

    int index = 0; 
    index = std::distance(offsets.begin(),
            std::find(offsets.begin(), offsets.end(), Point::Zeros()));
    UNIT_TEST((index < offsets.size()));
    UNIT_TEST((coefs[index] == 6.0));
    
    index = std::distance(offsets.begin(),
            std::find(offsets.begin(), offsets.end(), Point::Ones()));
    UNIT_TEST((index < offsets.size()));
    UNIT_TEST((coefs[index] == 2.0));
    
    index = std::distance(offsets.begin(),
            std::find(offsets.begin(), offsets.end(), Point::Basis(0,2)));
    UNIT_TEST((index < offsets.size()));
    UNIT_TEST((coefs[index] == 3.0));
    
    index = std::distance(offsets.begin(),
            std::find(offsets.begin(), offsets.end(), Point(6,5,4,3,2,1)));
    UNIT_TEST((index < offsets.size()));
    UNIT_TEST((coefs[index] == 4.0));
    END_TEST();
    
    //===================================
    BEGIN_TEST("Stencil Scalar Multiplication"); 
    Stencil<double> S0 = 1.0*Shift::Basis(0,-1)
                       - 2.0*Shift::Zeros()
                       + 1.0*Shift::Basis(0,1);
    auto S1 = 17.0*S0;
    if (VERBO > 1)
    {
        cout << "Initial Stencil S: " << endl;
        S0.print();
        cout << "S * 17.0: " << endl;
        S1.print();
    }
    
    std::vector<Point> offsets = S1.offsets(); 
    std::vector<double> coefs = S1.coefs();
    
    UNIT_TEST((offsets.size() == 3));
    UNIT_TEST((coefs.size() == 3));

    int index = 0; 
    index = std::distance(offsets.begin(),
            std::find(offsets.begin(), offsets.end(), Point::Basis(0,-1)));
    UNIT_TEST((index < offsets.size()));
    UNIT_TEST((coefs[index] == 17.0));
    
    index = std::distance(offsets.begin(),
            std::find(offsets.begin(), offsets.end(), Point::Zeros()));
    UNIT_TEST((index < offsets.size()));
    UNIT_TEST((coefs[index] == -34.0));
    
    index = std::distance(offsets.begin(),
            std::find(offsets.begin(), offsets.end(), Point::Basis(0,1)));
    UNIT_TEST((index < offsets.size()));
    UNIT_TEST((coefs[index] == 17.0));
     
    END_TEST();
    //===================================
    BEGIN_TEST("Stencil Composition"); 
    Stencil<double> S0 = 2.0*Shift::Basis(0);
    Stencil<double> S1 = 1.0*Shift::Basis(0,-1)
                       - 2.0*Shift::Zeros()
                       + 1.0*Shift::Basis(0,1);
    auto S2 = S0*S1;
    Stencil<double> S3 = 2.0*Shift::Zeros()
                       - 4.0*Shift::Basis(0,1)
                       + 2.0*Shift::Basis(0,2);
    if (VERBO > 1)
    {
        cout << "Composing the following two Stencils: " << endl;
        S0.print();
        S1.print();
        cout << "Resulting Stencil:" << endl;
        S2.print();
    }
    UNIT_TEST((S2 == S3));
    
    END_TEST();
    
    //===================================
    #if DIM > 1
    BEGIN_TEST("Transpose"); 
    Stencil<double> S0 = 1.0*Shift::Basis(1)
                       - 2.0*Shift::Zeros()
                       + 3.0*Shift::Basis(0,-1)
                       - 4.0*Shift::Ones();
    
    if (VERBO > 1)
    {
        cout << "Transposing the Stencil: " << endl;
        S0.print();
    }
    S0.transpose(0,1);
    if (VERBO > 1)
    {
        cout << "Resulting Stencil: " << endl;
        S0.print();
    }
    Stencil<double> S1 = 1.0*Shift::Basis(0)
                       - 2.0*Shift::Zeros()
                       + 3.0*Shift::Basis(1,-1)
                       - 4.0*Shift::Ones();
    UNIT_TEST((S0 == S1));
    END_TEST();
    #else
    OMIT_TEST("Stencil::transpose","Set DIM > 1 to run test");
    #endif
    
    //===================================
    BEGIN_TEST("Domain and Range"); 
    
    Bx B = Bx::Cube(16).shift(Point::Ones());
    Bx R, D;
    
    // 2*DIM + 1 Point Laplacian
    Stencil<double> S0 = (-2.0*DIM)*Shift::Zeros();
    for (int ii = 0; ii < DIM; ii++)
    {
        int d = ii+1;
        S0 += 1.0*Shift::Basis(ii,d);
        S0 += 1.0*Shift::Basis(ii,-d);
    }
    
    R = S0.range(B);
    D = S0.domain(B);
    if (VERBO > 1)
    {
        cout << "Computing Domain and Range of Stencil: " << endl;
        S0.print();
        cout << "Input: " << B << endl;
        cout << "Domain: " << D << endl;
        cout << "Range: " << R << endl;
    }
    UNIT_TEST((R == B.grow(Point(1,2,3,4,5,6)*-1)));
    UNIT_TEST((D == B.grow(Point(1,2,3,4,5,6))));

    //-----------------------------------
    // Linear Average
    Stencil<double> S1;
    Bx K = Bx::Cube(3);
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        S1 += 1.0*Shift(*iter);
    }
    S1.srcRatio() = Point::Ones(3);
    
    Bx r = Bx::Cube(3).shift(Point::Ones());
    Bx d = Bx(Point::Ones(2), Point::Ones(13));

    R = S1.range(d);
    D = S1.domain(r);
    if (VERBO > 1)
    {
        cout << "Computing Domain and Range of Stencil: " << endl;
        S1.print();
        cout << "Domain of " << r << ": " << D << endl;
        cout << "Range of " << d << ": " << R << endl;
    }
    UNIT_TEST((R == Bx(Point::Ones(), Point::Ones(3))));
    UNIT_TEST((D == Bx(Point::Ones(3), Point::Ones(11))));

    //-----------------------------------
    // Linear Interp
    Stencil<double> S2 = (2.0/3)*Shift::Zeros() + (1.0/3)*Shift::Basis(0);
    S2.destRatio() = Point::Ones(3);
    S2.destShift() = Point::Basis(0);
    
    r = Bx(Point::Ones(2), Point::Ones(12));
    d = Bx(Point::Ones(), Point::Ones(4));
    R = S2.range(d);
    D = S2.domain(r);
    if (VERBO > 1)
    {
        cout << "Computing Domain and Range of Stencil: " << endl;
        S2.print();
        cout << "Domain of " << r << ": " << D << endl;
        cout << "Range of " << d << ": " << R << endl;
    }
    UNIT_TEST((R == Bx(Point(4,3),Point(10,12))));    
    UNIT_TEST((D == Bx(Point::Ones(), Point::Ones(4))));
    END_TEST();

    //===================================
    BEGIN_TEST("Box Inference Logic");
    
    Point r_src = Point::Ones(2);
    Point r_dest = Point::Ones(3);

    Bx K = Bx::Cube(3).shift(Point::Ones(-1));
    Stencil<double> S;
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
      S += ((double)K.index(*iter))*Shift(*iter);
    }
    S.srcRatio() = r_src;
    S.destRatio() = r_dest;
    S.destShift() = Point::Basis(0);

    Bx B = Bx::Cube(7).shift(Point::Ones(-1));
    cout << B.taperCoarsen(r_src) << endl;

    END_TEST();

    

    //===================================
    BEGIN_TEST("Apply: Scalar Multiply");
    Stencil<double> S = 17.0*Shift::Zeros();
    Bx B = Bx::Cube(8);
    auto R = forall_p<double>([](Point p, Var<double> v)
    {
        v(0) = 1;
        for (int ii = 0; ii < DIM; ii++)
        {
            v(0) += p[ii];
        }
    },B);
    BoxData<double> D0 = S(R);
    Bx b = B.grow(-Point::Basis(0));
    BoxData<double> D1 = S(R,b);
    if (VERBO > 1)
    {
        #if DIM > 2
        cout << "Output is suppressed for DIM > 2" << endl;
        #else
        cout << "Applying Stencil: " << endl;
        S.print();
        cout << "Source Data: " << endl;
        R.print();
        cout << "Output with automatic Box: " << endl;
        D0.printData();
        cout << "Output on Bx " << b << endl;
        D1.printData();
        #endif
    }
    
    UNIT_TEST((D0.box() == B));
    UNIT_TEST((D1.box() == b));

    for (auto iter = B.begin(); iter != B.end(); ++iter)
    {
        UNIT_TEST((D0(*iter) == 17.0*R(*iter)));
        if (b.contains(*iter))
        {
            UNIT_TEST((D1(*iter) == 17.0*R(*iter)));
        }
    }
    END_TEST();
    
    //===================================
    BEGIN_TEST("Apply: Laplacian"); 
    Stencil<double> S0 = (-2.0*DIM)*Shift::Zeros();
    for (int ii = 0; ii < DIM; ii++)
    {
        S0 += 1.0*Shift::Basis(ii,+1);
        S0 += 1.0*Shift::Basis(ii,-1);
    }

    int numIter = 4;
    double error[numIter];
    
    int domainSize = 16;
    for (int ii = 0; ii < numIter; ii++)
    {
        double dx = M_PI/domainSize;
        Bx b = Bx::Cube(domainSize);
        Bx B = S0.domain(b);
        
        double ddx = dx*dx;
        auto S = S0*(1.0/ddx);

        auto R = forall_p<double>([dx](Point p, Var<double>& v)
        {
            v(0) = sin(p[0]*dx);
            #if DIM > 1
            v(0) += cos(p[1]*dx);
            #endif
        },B);
       
        BoxData<double> D0 = S(R);
        BoxData<double> D1 = S(R,b.grow(-Point::Basis(0)));
        BoxData<double> D2(B,1337);
        D2 |= S(R);
        BoxData<double> D3(B,17);
        D3 += S(R);
        
        UNIT_TEST((D0.box() == b));
        UNIT_TEST((D1.box() == b.grow(-Point::Basis(0))));
        UNIT_TEST((D2.box() == B));
        UNIT_TEST((D3.box() == B));
        
        for (auto iter = B.begin(); iter != B.end(); ++iter)
        {
            if (b.contains(*iter))
            {
                UNIT_TEST((abs(D2(*iter) - D0(*iter)) < 0.01));
                UNIT_TEST((abs(D3(*iter) - D0(*iter) - 17) < 0.01));
            } else {
                UNIT_TEST((D2(*iter) == 1337));
            }
            if (D1.box().contains(*iter))
            {
                UNIT_TEST((abs(D1(*iter) - D0(*iter)) < 0.01));
            }
        }

        D0 += R;
        error[ii] = D0.absMax();
        domainSize *= 2;
    }

    double rates[numIter-1];
    for (int ii = 1; ii < numIter; ii++)
    {
        rates[ii-1] = log2(error[ii-1]/error[ii]);
        if (VERBO > 1)
        {
            cout << "Error: " << error[ii] << " Rate: " << rates[ii-1] << endl; 
        }
    }
    UNIT_TEST((abs(rates[numIter-2] - 2) < 0.2));

    END_TEST();
    
    //===================================
    BEGIN_TEST("Apply: Average"); 
    Bx K = Bx::Cube(2);
    Stencil<double> S;
    double coef = 1.0/K.size();
    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        S += coef*Shift(*iter);
    }
    S.srcRatio() = Point::Ones(2);

    
    int domainSize = 16;
    Bx B0 = Bx::Cube(domainSize);
    Bx B1 = S.range(B0);

    auto Src = forall_p<double>([](Point p, Var<double>& v)
    {
        v(0) = 0;
        for (int ii = 0; ii < DIM; ii++)
        {
            v(0) += p[ii];
        }
    },B0);
   
    auto Soln = forall_p<double>([](Point p, Var<double>& v)
    {
        v(0) = 0;
        for (int ii = 0; ii < DIM; ii++)
        {
            v(0) += 2.0*p[ii];
        }
        v(0) += 1;
    },B1);
   
    BoxData<double> D0 = S(Src);
    BoxData<double> D1 = S(Src,B1.grow(-Point::Basis(0)));
    BoxData<double> D2(B1,1337);
    D2 |= S(Src);
    BoxData<double> D3(B1,17);
    D3 += S(Src);
    
    UNIT_TEST((D0.box() == B1));
    UNIT_TEST((D1.box() == B1.grow(-Point::Basis(0))));
    UNIT_TEST((D2.box() == B1));
    UNIT_TEST((D3.box() == B1));
    
    for (auto iter = B1.begin(); iter != B1.end(); ++iter)
    {
        UNIT_TEST((Soln(*iter) == D0(*iter)));
        UNIT_TEST((abs(D2(*iter) - D0(*iter)) < 1e-6));
        UNIT_TEST((abs(D3(*iter) - D0(*iter) - 17) < 1e-6));
        if (D1.box().contains(*iter))
        {
            UNIT_TEST((abs(D1(*iter) - D0(*iter)) < 1e-6));
        }
    }
   
    END_TEST();
    END_TEST_SUITE();
    
    //***********************************
    //  INTERPSTENCIL TESTS
    //***********************************
    BEGIN_TEST_SUITE(5,"Proto::InterpStencil");
    
    //===================================
    BEGIN_TEST("Default Constructor"); 
    InterpStencil<double> IS;
    UNIT_TEST((IS.empty()));
    UNIT_TEST((IS.kernel().size() <= 0));
    END_TEST();
    //===================================
    BEGIN_TEST("Standard Constructor");
    Point r = Point(2,3,4,5,6,7);
    InterpStencil<double> IS(r);
    UNIT_TEST((IS.kernel() == Bx(r - Point::Ones())));
    END_TEST();
    //===================================
    BEGIN_TEST("PiecewiseConstant and PiecewiseLinear");
    Point r = Point::Ones(3);
    auto PWC = InterpStencil<double>::PiecewiseConstant(r);
    auto PWL = InterpStencil<double>::PiecewiseLinear(r);
    
    int domainSize = 16;
    int numIter = 4;
    double error_C[numIter];
    double error_L[numIter];
    for (int ii = 0; ii < numIter; ii++)
    {
        Bx B0 = Bx::Cube(domainSize);
        Bx B1 = Bx(B0.low()*r, B0.high()*r);
        Bx B2 = B0.refine(r);
        BoxData<double> Src(B0);
        BoxData<double> DC0(B2,1337);
        BoxData<double> DL0(B1,1337);
        BoxData<double> DC1(B2,17);
        BoxData<double> DL1(B1,17);
        BoxData<double> Soln(B2);
        
        //double dx = (M_PI/4.0)/domainSize;
        double dx = 1.0/domainSize;
        forallInPlace_p([dx](Point p, Var<double>& v)
        {
            //v(0) = p[0]*dx;
            v(0) = sin(dx*p[0]);
            #if DIM > 1
            v(0) *= cos(dx*p[1]);
            #endif
            
        },Src);
        
        forallInPlace_p([dx](Point p, Var<double>& v)
        {
            //v(0) = p[0]*dx/3.0;
            v(0) = sin(dx/3.0*p[0]);
            #if DIM > 1
            v(0) *= cos(dx/3.0*p[1]);
            #endif
        },Soln);

        DC0 |= PWC(Src);
        DL0 |= PWL(Src);
        DC1 += PWC(Src);
        DL1 += PWL(Src);
        BoxData<double> DC2 = PWC(Src);
        BoxData<double> DL2 = PWL(Src);
    
        UNIT_TEST((DC2.box() == B2));
        UNIT_TEST((DL2.box() == B1));
        DC0 -= DC2;
        DL0 -= DL2;
        DC1 -= 17;
        DL1 -= 17;
        DC1 -= DC2;
        DL1 -= DL2;

        UNIT_TEST((DC0.absMax() < 1e-6));
        UNIT_TEST((DL0.absMax() < 1e-6));
        UNIT_TEST((DC1.absMax() < 1e-6));
        UNIT_TEST((DL1.absMax() < 1e-6));
        DC2 -= Soln;
        DL2 -= Soln;
        error_C[ii] = DC2.absMax();
        error_L[ii] = DL2.absMax();
        domainSize *= 2;
    }
    
    for (int ii = 1; ii < numIter; ii++)
    {
        double rate = log2(error_C[ii-1]/error_C[ii]);
        UNIT_TEST((abs(rate - 1) < 0.1));
        if (VERBO > 1)
        {
            cout << "PWC Error: " << error_C[ii] << " Rate: " << rate << endl;
        }
    }
    for (int ii = 1; ii < numIter; ii++)
    {
        double rate = log2(error_L[ii-1]/error_L[ii]);
        UNIT_TEST((abs(rate - 2) < 0.1));
        if (VERBO > 1)
        {
            cout << "PWL Error: " << error_L[ii] << " Rate: " << rate << endl;
        }
    }
    END_TEST();
    
    //===================================
    BEGIN_TEST("Box Inference");
    Bx B0 = Bx::Cube(4).grow(1);
    Bx B1 = Bx::Cube(8).grow(1);
    Bx K = Bx(Point::Ones(-2),Point::Ones(8));
    auto Src = forall_p<double>([](Point p, Var<double>& a_v)
    {
        a_v(0) = 0;
        for (int ii = 0; ii < DIM; ii++)
        {
            a_v(0) += p[ii];
        }
    },B0);
    auto Soln = forall_p<double>([](Point p, Var<double>& a_v)
    {
        a_v(0) = 0;
        for (int ii = 0; ii < DIM; ii++)
        {
            a_v(0) += p[ii]/2.0;
        }
    },K);
    
    BoxData<double> Dest0(B1,1337);
    BoxData<double> Dest1(B1,1337);
    BoxData<double> Dest2(B1,1337);
    auto I = InterpStencil<double>::PiecewiseLinear(Point::Ones(2));
    
    Dest0 |= I(Src, B0.grow(-1));
    Dest1 |= I(Src, B0.grow(1));
    Dest2 |= I(Src);
    BoxData<double> Dest3 = I(Src);
    if (VERBO > 1)
    {
        cout << "Source Data" << endl;
        Src.printData(1);
        cout << "Solution Data" << endl;
        Soln.printData(1);
        cout << "Interpolating with input Bx: " << B0.grow(-1) << endl;
        Dest0.printData(1);
        cout << "Interpolating with input Bx: " << B0.grow(1) << endl;
        Dest1.printData(1);
        cout << "Interpolating with automatic Bx" << endl;
        Dest2.printData(1);
        cout << "Interpolating with automatic Bx and Destination" << endl;
        Dest3.printData(1);
    }
   
    UNIT_TEST((Dest3.box() == K));

    for (auto iter = K.begin(); iter != K.end(); ++iter)
    {
        UNIT_TEST((Dest3(*iter) == Soln(*iter)));
        if (B1.grow(-1).contains(*iter))
        {
            UNIT_TEST((Dest0(*iter) == Soln(*iter)));
        }
        if (B1.contains(*iter))
        {
            UNIT_TEST((Dest1(*iter) == Soln(*iter)));
            UNIT_TEST((Dest2(*iter) == Soln(*iter)));
        }
    }
    
    END_TEST();
    
    //===================================
    BEGIN_TEST("Less Trivial Application");
    int domainSize = 8;
    
    int numIter = 5;
    for (int ii = 0; ii < numIter; ii++)
    {
        Bx B0 = Bx::Cube(domainSize).grow(1);
        Bx B1 = Bx::Cube(2*domainSize).grow(1);
        double dx = M_PI/domainSize;
        BoxData<double> Src = forall_p<double>([dx](Point p, Var<double>& a_v)
        {
            double x = p[0]*dx;
            double y = p[1]*dx;
            a_v(0) = cos(x)*cos(y);
        },B0);
        BoxData<double> Soln = forall_p<double>([dx](Point p, Var<double>& a_v)
        {
            double x = p[0]*dx/2.0;
            double y = p[1]*dx/2.0;
            a_v(0) = cos(x)*cos(y) + cos(x);
        },B1);
        BoxData<double> Dest = forall_p<double>([dx](Point p, Var<double>& a_v)
        {
            double x = p[0]*dx/2.0;
            a_v(0) = cos(x);
        },B1);
        auto interp = InterpStencil<double>::PiecewiseLinear(Point::Ones(2));
        Dest += interp(Src);
        for (auto iter = B1.begin(); iter != B1.end(); ++iter)
        {
            double x =(*iter)[0]*dx/2.0;
            double y =(*iter)[1]*dx/2.0;
            Point p = (*iter) % Point::Ones(2);
            Point q = (*iter) / Point::Ones(2);
            double value = cos(x);
            if ( p == Point::Zeros())
            {
                value += Src(q);
                UNIT_TEST((abs(Dest(*iter) - value) < 1e-15));
            }
            else if ((p == Point::Basis(0)) ||
                     (p == Point::Basis(1)) ||
                     (p == Point::Basis(0,-1)) ||
                     (p == Point::Basis(1,-1)))
                      
            {
                value += (Src(q) + Src(q+p))/2.0;
                UNIT_TEST((abs(Dest(*iter) - value) < 1e-15));
            }
            else if (p ==  Point::Ones())
            {
                value += (Src(q) + Src(q+p) + Src(q + Point::Basis(0)) + Src(q + Point::Basis(1)))/4.0;
                UNIT_TEST((abs(Dest(*iter) - value) < 1e-15));
            } 
        }
        Dest -= Soln;
        if (VERBO > 1)
        {
            cout << "Error: " << Dest.absMax() << endl;
        }
        domainSize *= 2;
    }
    END_TEST();

    END_TEST_SUITE();
    } //end loop over all tests
} // end main

