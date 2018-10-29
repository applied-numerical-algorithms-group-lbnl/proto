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
    a_didTestPass = true;
    a_errorCode   = 0;
  }
  /**/
}
