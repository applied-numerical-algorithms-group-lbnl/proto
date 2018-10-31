#include "UnitTestFunctions.H"
#include <iostream>
using std::cout;
using std::endl;
using namespace Proto;
typedef Var<double,DIM> V;


PROTO_KERNEL_START
void iotaFuncF(Point           & a_p,
               V               & a_X,
               double            a_h)
{
  for (int ii = 0; ii < DIM; ii++)
  {
    a_X(ii) = a_p[ii]*a_h;  //for some reason, this was written without the 0.5 orginally
  }
}
PROTO_KERNEL_END(iotaFuncF,iotaFunc)

PROTO_KERNEL_START
void fooFuncF(Var<double,DIM>&       a_v, 
              const V              & a_x, 
              const Var<int>&        a_c)
{
  for (int ii = 0; ii < DIM; ii++)
  {
    a_v(ii) = a_x(ii) + a_c(0);
  }
}
PROTO_KERNEL_END(fooFuncF,fooFunc)

PROTO_KERNEL_START
void squareFuncF(Point& a_p, Var<double>& a_v, const Var<double>& a_c)
{
  double x = a_p[0];
  a_v(0) = x*x + a_c(0);
}
PROTO_KERNEL_END(squareFuncF,squareFunc)

namespace prototest
{
/****/
  bool UNIT_TEST(bool a_test, int& a_errorCode, int a_testnumber)
  {
    if(!a_test)
    {
      a_errorCode = a_testnumber;
    }
    else
    {
      a_errorCode = 0;
    }
    return a_test;
  }
  void printTestMessage(const string& a_testName, const int & a_errorCode, bool a_didTestPass)
  {
    if(a_didTestPass)
    {
      cout << a_testName << " passed all tests" << endl;
    }
    else
    {
      cout << a_testName << " failed test with error code " << a_errorCode << endl;
    }
  }
/****/
  void pointTest(int& a_errorCode, bool & a_didTestPass)
  {
    //default construction
    {
      Point p;
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((p[ii]==0), a_errorCode, 1);
        if(!a_didTestPass) return;
      }
    }
    // array construction
    {
      int v[DIM];
      for (int ii = 0; ii < DIM; ii++)
      {
        v[ii] = (ii+1)*17;
      }
      Point p(v);
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((p[ii]==v[ii]), a_errorCode, 2);
        if(!a_didTestPass) return;
      }
    }
    //variadic construction
    {
      Point p(1,2,3,4,5,6);
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((p[ii]==(ii+1)), a_errorCode, 3);
        if(!a_didTestPass) return;
      }
    }
//copy construction
    {
      Point p(1,2,3,4,5,6);
      Point q(p);
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((p[ii]==q[ii]), a_errorCode, 4);
        if(!a_didTestPass) return;
        a_didTestPass = UNIT_TEST((&p[ii]!=&q[ii]), a_errorCode, 5);
        if(!a_didTestPass) return;
      }
    }
    //static methods
    {
      Point p0 = Point::Zeros();
      Point p1 = Point::Ones();
      Point p2 = Point::Ones(17);
      Point p3 = Point::Basis(0);
      Point p4 = Point::Basis(DIM-1,17);
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((p0[ii]==0), a_errorCode, 6);
        if(!a_didTestPass) return;
        a_didTestPass = UNIT_TEST((p1[ii]==1), a_errorCode, 7);
        if(!a_didTestPass) return;
        a_didTestPass = UNIT_TEST((p2[ii]==17), a_errorCode, 8);
        if(!a_didTestPass) return;
        if (ii == 0)
        {
          a_didTestPass = UNIT_TEST((p3[ii] == 1), a_errorCode, 9);
          if(!a_didTestPass) return;
        }
        else
        {
          a_didTestPass = UNIT_TEST((p3[ii] == 0), a_errorCode, 10);
          if(!a_didTestPass) return;
        }
        if (ii == DIM-1)
        {
          a_didTestPass = UNIT_TEST((p4[ii] == 17), a_errorCode, 11);
          if(!a_didTestPass) return;
        }
        else
        {
          a_didTestPass = UNIT_TEST((p4[ii] == 0), a_errorCode, 12);
          if(!a_didTestPass) return;
        }
      }
    }
    //Accessor Methods
    {
      //the 13 test was exactly the same as the 14
      const Point q(1,2,3,4,5,6);
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((q[ii] == ii+1), a_errorCode, 14);
        if(!a_didTestPass) return;
      }
    }
    //Algebraic Operators
    {
      Point p0(1,2,3,4,5,6);
      Point p1;
      int flip = -1;
      for (int ii = 0; ii < DIM; ii++)
      {
        p1[ii] = p0[ii]*17*flip;
        flip *= -1;
      }

      Point p2 = p0 + p1;
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((p2[ii] == (p0[ii] + p1[ii])), a_errorCode, 15);
        if(!a_didTestPass) return;
      }
    
      p2 = p1 - p0;
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((p2[ii] == (p1[ii] - p0[ii])), a_errorCode, 16);
        if(!a_didTestPass) return;
      }

      p2 = p1*p0;
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((p2[ii] == (p1[ii]*p0[ii])), a_errorCode, 17);
        if(!a_didTestPass) return;
      }

      p2 = p1/p0;
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((p2[ii] == (p1[ii]/p0[ii])), a_errorCode, 13); 
        if(!a_didTestPass) return;
      }

      p2 = p1/17;
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((p2[ii] == (p1[ii]/17)), a_errorCode, 18); 
        if(!a_didTestPass) return;
      }

      p1 = Point::Ones(2); 
      p2 = p0 % p1;
      for (int ii = 0; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((p2[ii] == (p0[ii] % p1[ii])), a_errorCode, 19); 
        if(!a_didTestPass) return;
      }
    }  
    a_didTestPass = true;
    a_errorCode   = 0;
  }
/****/
  void bxTest(int& a_errorCode, bool & a_didTestPass)
  {
    //Default Constructor
    {
      Bx B;
      a_didTestPass = UNIT_TEST((B.low() == Point::Zeros())  , a_errorCode, 20) ;
      if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((B.high() == Point::Ones(-1)), a_errorCode, 21) ;
      if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((B.size() <= 0)              , a_errorCode, 22) ;
      if(!a_didTestPass) return;
    }
    //Point Constructors
    {
      Point p1(1,2,3,4,5,6);
      Point p0 = (-1)*p1;
      Bx B0(p1,p0);
      Bx B1(p0,p1);
      Bx B2(p1);
      int s1 = 1; int s2 = 1;
      for (int ii = 0; ii < DIM; ii++)
      {
        s1 *= (p1[ii] - p0[ii] + 1);
        s2 *= (p1[ii] + 1);
      }
      a_didTestPass = UNIT_TEST((B0.size() <= 0)   , a_errorCode, 23) ;
      if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((B0.empty())       , a_errorCode, 24) ;
      if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((B1.size() == s1)  , a_errorCode, 25) ;
      if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((!B1.empty())      , a_errorCode, 26) ;
      if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((B2.size() == s2)  , a_errorCode, 27) ;
      if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((!B2.empty()    )  , a_errorCode, 28) ;
      if(!a_didTestPass) return;
    }
    //Copy Constructor
    {
      Bx B0(Point(1,2,3,4,5,6));
      Bx B1(B0);
      a_didTestPass = UNIT_TEST((B0==B1        ) , a_errorCode, 29) ;
      if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((&B0!=&B1      ) , a_errorCode, 30) ;
      if(!a_didTestPass) return;
    }
    //Static Methods
    {
      int size = 17;
      Bx B = Bx::Cube(size);
      a_didTestPass = UNIT_TEST((B.size() == ipow<DIM>(size))    , a_errorCode, 31); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((B.low() == Point::Zeros())      , a_errorCode, 32); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((B.high() == Point::Ones(size-1)), a_errorCode, 33); if(!a_didTestPass) return;
    }
    //Iteration and Indexing
    {
      Point max = Point(1,2,3,4,5,6);
      Point min = -1*max;
      Bx B(min,max);
      int index = 0;
      Point p = min;
      for (auto iter = B.begin(); iter != B.end(); ++iter, ++index)
      {
        a_didTestPass = UNIT_TEST((*iter == p)             , a_errorCode, 34); if(!a_didTestPass) return;
        a_didTestPass = UNIT_TEST((B.index(*iter) == index), a_errorCode, 35); if(!a_didTestPass) return;
        ++p[0];
        for (int ii = 0; ii < DIM-1; ii++)
        {
          if (p[ii] > max[ii])
          {
            p[ii] = min[ii];
            ++p[ii+1];
          }
        }
      }
      p = max;
      index = B.size()-1;
      for (auto iter = B.rbegin(); iter != B.rend(); --iter, --index)
      {
        a_didTestPass = UNIT_TEST((*iter == p)             , a_errorCode, 36); if(!a_didTestPass) return;
        a_didTestPass = UNIT_TEST((B.index(*iter) == index), a_errorCode, 37); if(!a_didTestPass) return;
        --p[0];
        for (int ii = 0; ii < DIM-1; ii++)
        {
          if (p[ii] < min[ii])
          {
            p[ii] = max[ii];
            --p[ii+1];
          }
        }
      }
    }
    //Shift
    {
      Bx b0 = Bx::Cube(16);
      Point s(1,-2,3,-4,5,-6);
    
      Bx b1 = b0.shift(s);

      a_didTestPass = UNIT_TEST(((b0.low() + s) == b1.low())  , a_errorCode, 38); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST(((b0.high() + s) == b1.high()), a_errorCode, 39); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((&b0 != &b1)                  , a_errorCode, 40); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b0 == Bx::Cube(16))          , a_errorCode, 41); if(!a_didTestPass) return;
    }
    //Grow
    {
      Bx b0 = Bx::Cube(16);
      Point s(1,-2,3,-4,5,-6);
      //Grow (Point)
      Bx b1 = b0.grow(s);

      a_didTestPass = UNIT_TEST(((b0.low() - s) == b1.low())  , a_errorCode, 42); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST(((b0.high() + s) == b1.high()), a_errorCode, 43); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((&b0 != &b1)                  , a_errorCode, 44); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b0 == Bx::Cube(16))          , a_errorCode, 45); if(!a_didTestPass) return;

      //Grow (scalar)
      b1 = b0.grow(3);

      a_didTestPass = UNIT_TEST(((b0.low() - 3) == b1.low())  , a_errorCode, 46); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST(((b0.high() + 3) == b1.high()), a_errorCode, 47); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((&b0 != &b1)                  , a_errorCode, 48); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b0 == Bx::Cube(16))          , a_errorCode, 49); if(!a_didTestPass) return;                               
    }
    //Coarsen
    {
      Point low = Point::Ones(-2);
      Point high = Point::Ones(3);
      Point r = Point::Ones() + Point::Basis(0); //(2,1,1,...,1)
      Bx b0 = Bx(low,high); 
      Bx b1 = b0.coarsen(2); 
      Bx b2 = b0.coarsen(3);
      Bx b3 = b0.coarsen(r);

      a_didTestPass = UNIT_TEST((b1.low() == b0.low()/2)          , a_errorCode, 50); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b1.high() == b0.high()/2)        , a_errorCode, 51); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b2.low() == b0.low()/3)          , a_errorCode, 52); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b2.high() == b0.high()/3)        , a_errorCode, 53); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b3.low()[0] == b0.low()[0]/2)    , a_errorCode, 54); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b3.high()[0] == b0.high()[0]/2)  , a_errorCode, 55); if(!a_didTestPass) return;

      for (int ii = 1; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((b3.low()[ii] == b0.low()[ii])  , a_errorCode, 56); if(!a_didTestPass) return;
        a_didTestPass = UNIT_TEST((b3.high()[ii] == b0.high()[ii]), a_errorCode, 57); if(!a_didTestPass) return;
      }

      a_didTestPass = UNIT_TEST((&b0 != &b1)        , a_errorCode, 58); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((&b0 != &b2)        , a_errorCode, 59); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((&b0 != &b3)        , a_errorCode, 60); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b0 == Bx(low,high)), a_errorCode, 61); if(!a_didTestPass) return;
    }
    //refine
    {
      Point low = Point::Ones(-2);
      Point high = Point::Ones(3);
      Point r = Point::Ones() + Point::Basis(0); //(2,1,1,...,1)
      Bx b0 = Bx(low,high);
      Bx b1 = b0.refine(2);
      Bx b2 = b0.refine(3);
      Bx b3 = b0.refine(r);
      a_didTestPass = UNIT_TEST((b1.low() == b0.low()*2)                                , a_errorCode, 62); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b1.high() == (b0.high()+Point::Ones())*2-Point::Ones()), a_errorCode, 63); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b2.low() == b0.low()*3)                                , a_errorCode, 64); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b2.high() == (b0.high()+Point::Ones())*3-Point::Ones()), a_errorCode, 65); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b3.low()[0] == b0.low()[0]*2)                          , a_errorCode, 66); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b3.high()[0] == (b0.high()[0]+1)*2 - 1)                , a_errorCode, 67); if(!a_didTestPass) return;

      for (int ii = 1; ii < DIM; ii++)
      {
        a_didTestPass = UNIT_TEST((b3.low()[ii] == b0.low()[ii])  , a_errorCode, 68); if(!a_didTestPass) return;
        a_didTestPass = UNIT_TEST((b3.high()[ii] == b0.high()[ii]), a_errorCode, 69); if(!a_didTestPass) return;
      }
      a_didTestPass = UNIT_TEST((&b0 != &b1)        , a_errorCode, 70); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((&b0 != &b2)        , a_errorCode, 71); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((&b0 != &b3)        , a_errorCode, 72); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((b0 == Bx(low,high)), a_errorCode, 73); if(!a_didTestPass) return;

    }
    //flatten
    {
      Bx b0 = Bx::Cube(17);
      for (int ii = 0; ii < DIM; ii++)
      {
        Bx b1 = b0.flatten(ii);
        Bx b2 = b0.flatten(ii,true);
        a_didTestPass = UNIT_TEST((b1.low()  == b0.low()) , a_errorCode, 74); if(!a_didTestPass) return;
        a_didTestPass = UNIT_TEST((b2.high() == b0.high()), a_errorCode, 75); if(!a_didTestPass) return;
        for (int jj = 0; jj < DIM; jj++)
        {
          if (jj == ii)
          {
            a_didTestPass = UNIT_TEST((b1.high()[jj] == b1.low()[jj]) , a_errorCode, 76); if(!a_didTestPass) return;
            a_didTestPass = UNIT_TEST((b2.low()[jj]  == b2.high()[jj]), a_errorCode, 77); if(!a_didTestPass) return;
          }   
          else
          {   
            a_didTestPass = UNIT_TEST((b1.high()[jj] == b0.high()[jj]), a_errorCode, 78); if(!a_didTestPass) return;
            a_didTestPass = UNIT_TEST((b2.low()[jj]  == b0.low()[jj]) , a_errorCode, 79); if(!a_didTestPass) return;
          }
        }
      }
    }
    //extrude
    {
      Bx b0 = Bx::Cube(17);
      for (int ii = 0; ii < DIM; ii++)
      {
        Bx b1 = b0.extrude(ii,3);
        Bx b2 = b0.extrude(ii,3,true);
        a_didTestPass = UNIT_TEST((b1.high() == b0.high()), a_errorCode, 80); if(!a_didTestPass) return;
        a_didTestPass = UNIT_TEST((b2.low()  == b0.low()) , a_errorCode, 81); if(!a_didTestPass) return;
        for (int jj = 0; jj < DIM; jj++)
        {
          if (jj == ii)
          {
            a_didTestPass = UNIT_TEST((b1.low()[jj]  == b0.low()[jj]-3) , a_errorCode, 82); if(!a_didTestPass) return;
            a_didTestPass = UNIT_TEST((b2.high()[jj] == b0.high()[jj]+3), a_errorCode, 83); if(!a_didTestPass) return;
          }                                                               
          else                                                            
          {                                                               
            a_didTestPass = UNIT_TEST((b1.low()[jj]  == b0.low()[jj])   , a_errorCode, 84); if(!a_didTestPass) return;
            a_didTestPass = UNIT_TEST((b2.high()[jj] == b0.high()[jj])  , a_errorCode, 85); if(!a_didTestPass) return;
          }
        }
      }
    }
    //Mod
    {
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
    
      a_didTestPass = UNIT_TEST((q0 == p0)  , a_errorCode, 86); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((q1 == p1)  , a_errorCode, 87); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((q2 == low) , a_errorCode, 88); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((q3 == high), a_errorCode, 89); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((q4 == high), a_errorCode, 90); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((q5 == low) , a_errorCode, 91); if(!a_didTestPass) return;
    
    }

    a_didTestPass = true;
    a_errorCode   = 0;
  }
  /**/
  void boxdataTest(int& a_errorCode, bool & a_didTestPass)
  {
    //default constructor
    {
      BoxData<double,2,3> BD;
    
      a_didTestPass = UNIT_TEST((BD.box() == Bx(Point::Ones(-1))), a_errorCode, 100); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((BD.size() <= 0                 ), a_errorCode, 101); if(!a_didTestPass) return;
    }
    //Box Constructor
    {
      Bx B = Bx(Point(1,2,3,4,5,6,7));
      BoxData<int, 3, 4, 5> BD(B);
      int size = 1; 
      for (int ii = 0; ii < DIM; ii++)
      {
        size *= (ii+2);
      }
      a_didTestPass = UNIT_TEST((BD.box() == B          ), a_errorCode, 102); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((BD.size() == size*3*4*5), a_errorCode, 103); if(!a_didTestPass) return;
    }
    //("Initialization Constructor (And Accessor)"); 
    {
      Bx B = Bx(Point(1,2,3,4,5,6,7));
      BoxData<int, 3, 4, 5> BD(B,1337);
      int size = 1; 
      for (int ii = 0; ii < DIM; ii++)
      {
        size *= (ii+2);
      }
      a_didTestPass = UNIT_TEST((BD.box() == B)          , a_errorCode, 104); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((BD.size() == size*3*4*5), a_errorCode, 105); if(!a_didTestPass) return;
      int maxval = BD.max();      
      int minval = BD.min();
      a_didTestPass = UNIT_TEST((maxval == 1337), a_errorCode, 106); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((minval == 1337), a_errorCode, 107); if(!a_didTestPass) return;

    }
    //Move constructor
    {
#ifndef PROTO_MEM_CHECK
      cout << "Move Constructor test omittedCompile with PROTO_MEM_CHECK=TRUE to run this test" << endl;
#else
      memcheck::FLUSH_CPY();
      Bx B = Bx(Point(1,2,3,4,5,6,7)*2);
      double dx = 0.1;
      auto X = forall_p<double,DIM>(iotaFunc,B,dx);
//    int ncpy1 = memcheck::numcopies;
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0),  a_errorCode, 108); if(!a_didTestPass) return;
      BoxData<double,DIM> Y(B,1337);
      memcheck::FLUSH_CPY();
      Y = forall_p<double,DIM>(iotaFunc,B,dx);
//    int ncpy2 = memcheck::numcopies;
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 109); if(!a_didTestPass) return;

      BoxData<int> errf = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM> xv,
                                                          Var<double, DIM> yv) 
                                        {  
                                          err(0) = 0;
                                          for (int ii = 0; ii < DIM; ii++)
                                          {
                                            if(xv(ii) != dx*p[ii])
                                            {
                                              err(0) = 1;
                                            }
                                            if(yv(ii) != dx*p[ii])
                                            {
                                              err(0) = 2;
                                            }
                                          }
                                        }, B, X, Y);
      a_didTestPass = UNIT_TEST((errf.max() == 0), a_errorCode, 110); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errf.min() == 0), a_errorCode, 111); if(!a_didTestPass) return;

#endif
    }
    //Cinterval
    {

      CInterval I0(1,2,3,4,5,6);
      CInterval I1{{1,2},{3,4},{5,6}};
      CInterval I2{{},{3,4},{}};
      CInterval I3{1,2};

      a_didTestPass = UNIT_TEST((I0.low(0) == 1 ), a_errorCode, 112); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I0.high(0) == 2), a_errorCode, 113); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I0.low(1) == 3 ), a_errorCode, 114); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I0.high(1) == 4), a_errorCode, 115); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I0.low(2) == 5 ), a_errorCode, 116); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I0.high(2) == 6), a_errorCode, 117); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I1.low(0) == 1 ), a_errorCode, 118); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I1.high(0) == 2), a_errorCode, 119); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I1.low(1) == 3 ), a_errorCode, 120); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I1.high(1) == 4), a_errorCode, 121); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I1.low(2) == 5 ), a_errorCode, 122); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I1.high(2) == 6), a_errorCode, 123); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I2.low(0) == 0 ), a_errorCode, 124); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I2.high(0) == 0), a_errorCode, 125); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I2.low(1) == 3 ), a_errorCode, 126); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I2.high(1) == 4), a_errorCode, 127); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I2.low(2) == 0 ), a_errorCode, 128); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I2.high(2) == 0), a_errorCode, 129); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I3.low(0) == 1 ), a_errorCode, 130); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I3.high(0) == 2), a_errorCode, 131); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I3.low(1) == 0 ), a_errorCode, 132); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I3.high(1) == 0), a_errorCode, 133); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I3.low(2) == 0 ), a_errorCode, 134); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((I3.high(2) == 0), a_errorCode, 135); if(!a_didTestPass) return;

    }
    //copy
    {
#ifndef PROTO_MEM_CHECK
      cout << "Copy test omitted.  Compile with PROTO_MEM_CHECK=TRUE to run this test" << endl;
#else
      Bx B = Bx(Point(1,2,3,4,5,6,7)*2);
      Point s = Point::Ones();
      Bx B1 = B.shift(s);
      Bx b = B.grow(-1);
      double dx = 0.1;
      auto X = forall_p<double,DIM>(iotaFunc, B, dx);
    
      BoxData<double,DIM> Y0(B,1337);
      BoxData<double,DIM> Y1(B1,1337);

      memcheck::FLUSH_CPY();

      X.copyTo(Y0);

      a_didTestPass = UNIT_TEST((memcheck::numcopies == 1), a_errorCode, 136); if(!a_didTestPass) return;

      BoxData<int> errf = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM> xv,
                                                          Var<double, DIM> yv) 
                                        {  
                                          err(0) = 0;
                                          for (int ii = 0; ii < DIM; ii++)
                                          {
                                            if(xv(ii) != yv(ii))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, B, X, Y0);
      a_didTestPass = UNIT_TEST((errf.max() == 0), a_errorCode, 137); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errf.min() == 0), a_errorCode, 138); if(!a_didTestPass) return;

    
      memcheck::FLUSH_CPY();

      X.copyTo(Y1);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 1), a_errorCode, 139); if(!a_didTestPass) return;

      Bx Binter = B & B1;
      BoxData<int> errg = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM> xv,
                                                          Var<double, DIM> yv) 
                                        {  
                                          err(0) = 0;
                                          for (int ii = 0; ii < DIM; ii++)
                                          {
                                            if(xv(ii) != yv(ii))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, Binter, X, Y1);
      a_didTestPass = UNIT_TEST((errg.max() == 0), a_errorCode, 140); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errg.min() == 0), a_errorCode, 141); if(!a_didTestPass) return;
//left out the wacky copy-shift test as I could not quite figure out what it was doing.
#endif
    }
    //algebraic operations
    {
      Bx B0 = Bx::Cube(4);
//    Bx B1 = B0.shift(Point::Ones());
      Bx B1 = B0;
      double dx = 0.1;
      auto D0 = forall_p<double,DIM>(iotaFunc,B0,dx);
      BoxData<double,DIM> delta(B0,dx/2);
      D0 += delta;
      BoxData<double,DIM> D1(B1,17);

      D1 += D0;

      BoxData<int> errh = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM> xv,
                                                          Var<double, DIM> yv) 
                                        {  
                                          err(0) = 0;
                                          for (int ii = 0; ii < DIM; ii++)
                                          {
                                            if(xv(ii) != (yv(ii)+17))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, B0, D1, D0);
      a_didTestPass = UNIT_TEST((errh.max() == 0), a_errorCode, 142); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errh.min() == 0), a_errorCode, 143); if(!a_didTestPass) return;

      D1.setVal(17);
      D1 -= D0;


      BoxData<int> erri = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM> xv,
                                                          Var<double, DIM> yv) 
                                        {  
                                          err(0) = 0;
                                          for (int ii = 0; ii < DIM; ii++)
                                          {
                                            if(xv(ii) != (17-yv(ii)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, B0, D1, D0);
      a_didTestPass = UNIT_TEST((erri.max() == 0), a_errorCode, 144); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((erri.min() == 0), a_errorCode, 145); if(!a_didTestPass) return;

      D1.setVal(17);

      D1 *= D0;
      BoxData<int> errj = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM> xv,
                                                          Var<double, DIM> yv) 
                                        {  
                                          err(0) = 0;
                                          for (int ii = 0; ii < DIM; ii++)
                                          {
                                            if(xv(ii) != (17*yv(ii)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, B0, D1, D0);
      a_didTestPass = UNIT_TEST((errj.max() == 0), a_errorCode, 146); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errj.min() == 0), a_errorCode, 147); if(!a_didTestPass) return;

      D1.setVal(17);
      D1 /= D0;

      BoxData<int> errk = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM> xv,
                                                          Var<double, DIM> yv) 
                                        {  
                                          err(0) = 0;
                                          for (int ii = 0; ii < DIM; ii++)
                                          {
                                            if(xv(ii) != (17/yv(ii)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, B0, D1, D0);
      a_didTestPass = UNIT_TEST((errk.max() == 0), a_errorCode, 148); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errk.min() == 0), a_errorCode, 149); if(!a_didTestPass) return;
    }

    //reductions 
    {
      Bx B = Bx::Cube(4).shift(Point::Basis(0,-2));
      double dx = 1;
      auto D = forall_p<double,DIM>(iotaFunc,B,dx);

      a_didTestPass = UNIT_TEST((D.max()    ==  3), a_errorCode, 150); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((D.min()    == -2), a_errorCode, 151); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((D.absMax() ==  3), a_errorCode, 152); if(!a_didTestPass) return;
   
      for (int ii = 0; ii < DIM; ii++)
      {
        if (ii == 0)
        {
          a_didTestPass = UNIT_TEST((D.max(ii)    ==  1), a_errorCode, 153); if(!a_didTestPass) return;
          a_didTestPass = UNIT_TEST((D.min(ii)    == -2), a_errorCode, 154); if(!a_didTestPass) return;
          a_didTestPass = UNIT_TEST((D.absMax(ii) ==  2), a_errorCode, 155); if(!a_didTestPass) return;
        } 
        else 
        {
          a_didTestPass = UNIT_TEST((D.max(ii)    == 3), a_errorCode, 156); if(!a_didTestPass) return;
          a_didTestPass = UNIT_TEST((D.min(ii)    == 0), a_errorCode, 157); if(!a_didTestPass) return;
          a_didTestPass = UNIT_TEST((D.absMax(ii) == 3), a_errorCode, 158); if(!a_didTestPass) return;
        }
      }
    }
    //alias and slice
    {
#ifndef PROTO_MEM_CHECK
      cout << " omitting Alias and Slice test. To run this test, compile with PROTO_MEM_CHECK=TRUE" << endl;
#else
      Bx B0 = Bx::Cube(4).shift(Point::Basis(0,-2));
      Point shift = Point::Basis(0,-1);
      double dx = 0.1;
      auto D0 = forall_p<double,DIM>(iotaFunc,B0,dx);
    
      memcheck::FLUSH_CPY();
      auto D1 = alias(D0,shift);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 159); if(!a_didTestPass) return;
      
      //again, not so sure how to write this with the shifts
//    
//      for (int ii = 0; ii < DIM; ii++)
//        for (auto iter = B0.begin(); iter != B0.end(); ++iter)
//        {
//          a_didTestPass = UNIT_TEST((D1(*iter+shift,ii) == D0(*iter,ii)));
//          UNIT_TEST((&D1(*iter+shift,ii) == &D0(*iter,ii)));
//        }

      memcheck::FLUSH_CPY();
      auto D2 = slice(D0,1);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 160); if(!a_didTestPass) return;

    
      BoxData<int> errl = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, 1> xv,
                                                          Var<double, DIM> yv) 
                                        {  
                                          err(0) = 0;
                                          if(xv(0) != (yv(1)))
                                          {
                                            err(0) = 1;
                                          }
                                        }, B0, D2, D0);
      a_didTestPass = UNIT_TEST((errl.max() == 0), a_errorCode, 161); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errl.min() == 0), a_errorCode, 162); if(!a_didTestPass) return;

#endif
    }
    //=====================START FORALL==============
    {
#ifndef PROTO_MEM_CHECK
      cout << " omitting forall  test. To run this test, compile with PROTO_MEM_CHECK=TRUE" << endl;
#else
    
      const Bx B0 = Bx::Cube(5).shift(Point::Basis(0,-2));
      const Bx B1 = Bx::Cube(5);
      const Bx B2 = B0 & B1;
      const Bx b2 = Bx::Cube(2);
      double dx = 0.1;
      auto X = forall_p<double,DIM>(iotaFunc,B0,dx);
      BoxData<int> C(B1,17);

      // forall
      //-------------------------------------------
      // with automatic Box

      memcheck::FLUSH_CPY();
      Bx interbox = X.box() & C.box();
      BoxData<double,DIM> D0 = forall<double,DIM>(fooFunc,X,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 163); if(!a_didTestPass) return;
    
      a_didTestPass = UNIT_TEST((D0.box() == B2), a_errorCode, 164); if(!a_didTestPass) return;

      BoxData<int> errm = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM>    dv,
                                                          Var<double, DIM>    xv,
                                                          Var<int,      1>    cv) 
                                        {  
                                          err(0) = 0;
                                          for (int dir = 0; dir < DIM; dir++)
                                          {
                                            if(dv(dir) != (cv(0) + xv(dir)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, B2, D0, X, C);
      a_didTestPass = UNIT_TEST((errm.max() == 0), a_errorCode, 165); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errm.min() == 0), a_errorCode, 166); if(!a_didTestPass) return;
    
      // with supplied Box

      memcheck::FLUSH_CPY();
      BoxData<double,DIM> D1 = forall<double,DIM>(fooFunc,b2,X,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 167); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((D1.box() == b2), a_errorCode, 168); if(!a_didTestPass) return;

      BoxData<int> errn = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM>    dv,
                                                          Var<double, DIM>    xv,
                                                          Var<int,      1>    cv) 
                                        {  
                                          err(0) = 0;
                                          for (int dir = 0; dir < DIM; dir++)
                                          {
                                            if(dv(dir) != (cv(0) + xv(dir)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, b2, D1, X, C);
      a_didTestPass = UNIT_TEST((errn.max() == 0), a_errorCode, 169); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errn.min() == 0), a_errorCode, 170); if(!a_didTestPass) return;

      //forallInPlace
      //-------------------------------------------
      // with automatic box
      BoxData<double,DIM> D2(B1,1337);
      BoxData<double,DIM> D3(B1,1337);

      memcheck::FLUSH_CPY();
      forallInPlace(fooFunc,D2,X,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 171); if(!a_didTestPass) return;

      BoxData<int> erro = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM>    dv,
                                                          Var<double, DIM>    xv,
                                                          Var<int,      1>    cv) 
                                        {  
                                          err(0) = 0;
                                          for (int dir = 0; dir < DIM; dir++)
                                          {
                                            if(dv(dir) != (cv(0) + xv(dir)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, B2, D2, X, C);
      a_didTestPass = UNIT_TEST((erro.max() == 0), a_errorCode, 172); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((erro.min() == 0), a_errorCode, 173); if(!a_didTestPass) return;

      memcheck::FLUSH_CPY();
      forallInPlace(fooFunc,b2,D3,X,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 174); if(!a_didTestPass) return;
    
      BoxData<int> errp = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM>    dv,
                                                          Var<double, DIM>    xv,
                                                          Var<int,      1>    cv) 
                                        {  
                                          err(0) = 0;
                                          for (int dir = 0; dir < DIM; dir++)
                                          {
                                            if(dv(dir) != (cv(0) + xv(dir)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, b2, D3, X, C);
      a_didTestPass = UNIT_TEST((errp.max() == 0), a_errorCode, 175); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errp.min() == 0), a_errorCode, 176); if(!a_didTestPass) return;

      memcheck::FLUSH_CPY();
      forallInPlace(fooFunc,b2,D3,X,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 177); if(!a_didTestPass) return;
      BoxData<int> errq = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM>    dv,
                                                          Var<double, DIM>    xv,
                                                          Var<int,      1>    cv) 
                                        {  
                                          err(0) = 0;
                                          for (int dir = 0; dir < DIM; dir++)
                                          {
                                            if(dv(dir) != (cv(0) + xv(dir)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, b2, D3, X, C);

      a_didTestPass = UNIT_TEST((errq.max() == 0), a_errorCode, 177); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errq.min() == 0), a_errorCode, 178); if(!a_didTestPass) return;
#endif
    }    
    //forall_p==================================
    {
#ifndef PROTO_MEM_CHECK
      cout << "Forall_p test omitted. To run this test, please compile with PROTO_MEM_CHECK=TRUE" << endl;
#else
    
      Bx B0 = Bx::Cube(8);
      Bx B1 = Bx::Cube(8).shift(Point::Basis(0,-1));
      Bx B2 = B0 & B1;
      Bx b2 = B2.grow(-1);
    
      BoxData<double> C(B0,0.17);
    
      memcheck::FLUSH_CPY();
      auto D0 = forall_p<double>(squareFunc,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 179); if(!a_didTestPass) return;
    
      a_didTestPass = UNIT_TEST((D0.box() == B0), a_errorCode, 180); if(!a_didTestPass) return;
    
      BoxData<int> errr = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, 1>    dv)
                                        {  
                                          err(0) = 0;
                                          if(dv(0) != (p[0]*p[0]+0.17))
                                          {
                                            err(0) = 1;
                                          }
                                        }, B0, D0);
      a_didTestPass = UNIT_TEST((errr.max() == 0), a_errorCode, 179); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errr.min() == 0), a_errorCode, 180); if(!a_didTestPass) return;
      memcheck::FLUSH_CPY();
      auto D3 = forall_p<double>(squareFunc,b2,C);

      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 181); if(!a_didTestPass) return;
    
      a_didTestPass = UNIT_TEST((D3.box() == b2), a_errorCode, 182); if(!a_didTestPass) return;
      BoxData<int> errs = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, 1>    dv)
                                        {  
                                          err(0) = 0;
                                          if(dv(0) != (p[0]*p[0]+0.17))
                                          {
                                            err(0) = 1;
                                          }
                                        }, b2, D3);
      a_didTestPass = UNIT_TEST((errs.max() == 0), a_errorCode, 183); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errs.min() == 0), a_errorCode, 184); if(!a_didTestPass) return;
    

      BoxData<double> D1(B1,1337);
      BoxData<double> D2(B1,1337);

      memcheck::FLUSH_CPY();
      forallInPlace_p(squareFunc,D1,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 185); if(!a_didTestPass) return;
      BoxData<int> errt = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, 1>    dv)
                                        {  
                                          err(0) = 0;
                                          if(dv(0) != (p[0]*p[0]+0.17))
                                          {
                                            err(0) = 1;
                                          }
                                        }, B1 & B2, D1);
      a_didTestPass = UNIT_TEST((errt.max() == 0), a_errorCode, 186); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errt.min() == 0), a_errorCode, 187); if(!a_didTestPass) return;

    
      memcheck::FLUSH_CPY();
      forallInPlace_p(squareFunc,b2,D2,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 188); if(!a_didTestPass) return;
      BoxData<int> erru = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, 1>    dv)
                                        {  
                                          err(0) = 0;
                                          if(dv(0) != (p[0]*p[0]+0.17))
                                          {
                                            err(0) = 1;
                                          }
                                        }, B1 & b2, D2);
      a_didTestPass = UNIT_TEST((erru.max() == 0), a_errorCode, 189); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((erru.min() == 0), a_errorCode, 190); if(!a_didTestPass) return;
#endif
     
    }

    //forall random
    {
#ifndef PROTO_MEM_CHECK
      cout << "Forall box test omitted.  To run this test, please compile with PROTO_MEM_CHECK=TRUE " << endl;
#else
      const Bx B0 = Bx::Cube(5).shift(Point::Basis(0,-2));
      const Bx B1 = Bx::Cube(5);
      const Bx B2 = B0 & B1;
      const Bx b2 = Bx::Cube(2);
      double dx = 0.1;
      auto X = forall_p<double,DIM>(iotaFunc,B0,dx);
      BoxData<int> C(B1,17);
      // forall with automatic Box


      memcheck::FLUSH_CPY();
      BoxData<double,DIM> D0 = forall<double,DIM>(fooFunc,X,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 191); if(!a_didTestPass) return;
    
      a_didTestPass = UNIT_TEST((D0.box() == B2), a_errorCode, 192); if(!a_didTestPass) return;

      BoxData<int> errr = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM>    dv,
                                                          Var<double, DIM>    xv,
                                                          Var<int,      1>    cv) 
                                        {  
                                          err(0) = 0;
                                          for (int dir = 0; dir < DIM; dir++)
                                          {
                                            if(dv(dir) != (cv(0) + xv(dir)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, B2, D0, X, C);

      a_didTestPass = UNIT_TEST((errr.max() == 0), a_errorCode, 193); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errr.min() == 0), a_errorCode, 194); if(!a_didTestPass) return;
    
      // with supplied Box

      memcheck::FLUSH_CPY();
      BoxData<double,DIM> D1 = forall<double,DIM>(fooFunc,b2,X,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 195); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((D1.box() == b2), a_errorCode, 196); if(!a_didTestPass) return;
    
      BoxData<int> errs = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM>    dv,
                                                          Var<double, DIM>    xv,
                                                          Var<int,      1>    cv) 
                                        {  
                                          err(0) = 0;
                                          for (int dir = 0; dir < DIM; dir++)
                                          {
                                            if(dv(dir) != (cv(0) + xv(dir)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, b2, D1, X, C);

      a_didTestPass = UNIT_TEST((errs.max() == 0), a_errorCode, 197); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errs.min() == 0), a_errorCode, 198); if(!a_didTestPass) return;


      //forallInPlace
      //-------------------------------------------

      // with automatic box
      BoxData<double,DIM> D2(B1,1337);
      BoxData<double,DIM> D3(B1,1337);

      memcheck::FLUSH_CPY();
      forallInPlace(fooFunc,D2,X,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 199); if(!a_didTestPass) return;


      BoxData<int> errt = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM>    dv,
                                                          Var<double, DIM>    xv,
                                                          Var<int,      1>    cv) 
                                        {  
                                          err(0) = 0;
                                          for (int dir = 0; dir < DIM; dir++)
                                          {
                                            if(dv(dir) != (cv(0) + xv(dir)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, B1 & B2, D2, X, C);
      a_didTestPass = UNIT_TEST((errt.max() == 0), a_errorCode, 200); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errt.min() == 0), a_errorCode, 201); if(!a_didTestPass) return;

      memcheck::FLUSH_CPY();
      forallInPlace(fooFunc,b2,D3,X,C);

      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 202); if(!a_didTestPass) return;
      BoxData<int> erru = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, DIM>    dv,
                                                          Var<double, DIM>    xv,
                                                          Var<int,      1>    cv) 
                                        {  
                                          err(0) = 0;
                                          for (int dir = 0; dir < DIM; dir++)
                                          {
                                            if(dv(dir) != (cv(0) + xv(dir)))
                                            {
                                              err(0) = 1;
                                            }
                                          }
                                        }, B1 & b2, D3, X, C);
      a_didTestPass = UNIT_TEST((erru.max() == 0), a_errorCode, 200); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((erru.min() == 0), a_errorCode, 201); if(!a_didTestPass) return;

#endif
    }
    //==forall_p i guess=================================
    {    
#ifndef PROTO_MEM_CHECK
      cout << "omitting test Forall_p, To run this test, please compile with PROTO_MEM_CHECK=TRUE" << endl;
#else
      Bx B0 = Bx::Cube(8);
      Bx B1 = Bx::Cube(8).shift(Point::Basis(0,-1));
      Bx B2 = B0 & B1;
      Bx b2 = B2.grow(-1);
    
      BoxData<double> C(B0,0.17);
    
    
      memcheck::FLUSH_CPY();
      auto D0 = forall_p<double>(squareFunc,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 202); if(!a_didTestPass) return;

      a_didTestPass = UNIT_TEST((D0.box() == B0), a_errorCode, 203); if(!a_didTestPass) return;
    
      BoxData<int> errv = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, 1>    dv)
                                        {  
                                          err(0) = 0;
                                          if(dv(0) != (p[0]*p[0]+0.17))
                                          {
                                            err(0) = 1;
                                          }
                                        }, B0, D0);
      a_didTestPass = UNIT_TEST((errv.max() == 0), a_errorCode, 204); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errv.min() == 0), a_errorCode, 205); if(!a_didTestPass) return;
    
      memcheck::FLUSH_CPY();
      auto D3 = forall_p<double>(squareFunc,b2,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 206); if(!a_didTestPass) return;
    
      a_didTestPass = UNIT_TEST((D3.box() == b2), a_errorCode, 207); if(!a_didTestPass) return;
    
      BoxData<int> errw = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, 1>    dv)
                                        {  
                                          err(0) = 0;
                                          if(dv(0) != (p[0]*p[0]+0.17))
                                          {
                                            err(0) = 1;
                                          }
                                        }, b2, D3);
      a_didTestPass = UNIT_TEST((errw.max() == 0), a_errorCode, 208); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errw.min() == 0), a_errorCode, 209); if(!a_didTestPass) return;


      BoxData<double> D1(B1,1337);
      BoxData<double> D2(B1,1337);

      memcheck::FLUSH_CPY();
      forallInPlace_p(squareFunc,D1,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 210); if(!a_didTestPass) return;

      BoxData<int> errx = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, 1>    dv)
                                        {  
                                          err(0) = 0;
                                          if(dv(0) != (p[0]*p[0]+0.17))
                                          {
                                            err(0) = 1;
                                          }
                                        }, B1 & B2, D1);
      a_didTestPass = UNIT_TEST((errx.max() == 0), a_errorCode, 211); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errx.min() == 0), a_errorCode, 212); if(!a_didTestPass) return;
    
      memcheck::FLUSH_CPY();
      forallInPlace_p(squareFunc,b2,D2,C);
      a_didTestPass = UNIT_TEST((memcheck::numcopies == 0), a_errorCode, 213); if(!a_didTestPass) return;

      BoxData<int> errz = forall_p<int>([=] PROTO_LAMBDA (Point p, Var<int, 1> err, 
                                                          Var<double, 1>    dv)
                                        {  
                                          err(0) = 0;
                                          if(dv(0) != (p[0]*p[0]+0.17))
                                          {
                                            err(0) = 1;
                                          }
                                        }, B1 & b2, D2);
      a_didTestPass = UNIT_TEST((errz.max() == 0), a_errorCode, 213); if(!a_didTestPass) return;
      a_didTestPass = UNIT_TEST((errz.min() == 0), a_errorCode, 214); if(!a_didTestPass) return;
    

#endif
    }     
            
    a_didTestPass = true;
    a_errorCode   = 0;
  }

  /**/
}
