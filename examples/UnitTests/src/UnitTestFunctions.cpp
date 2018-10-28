#include "UnitTestFunctions.H"
#include <iostream>
using std::cout;
using std::endl;
namespace prototest
{
  using namespace Proto;
/****/
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
        if(p[ii] != 0)
        {
          a_didTestPass = false;
          a_errorCode = 1;
          return;
        }
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
        if(p[ii] != v[ii])
        {
          a_didTestPass = false;
          a_errorCode = 2;
          return;
        }
      }
    }
    //variadic construction
    {
      Point p(1,2,3,4,5,6);
      for (int ii = 0; ii < DIM; ii++)
      {
        if(p[ii] != ii+1)
        {
          a_didTestPass = false;
          a_errorCode = 3;
          return;
        }
      }
    
      a_didTestPass = true;
      a_errorCode   = 0;
    }
//copy construction
    {
      Point p(1,2,3,4,5,6);
      Point q(p);
      for (int ii = 0; ii < DIM; ii++)
      {
        if(p[ii] != q[ii])
        {
          a_didTestPass = false;
          a_errorCode = 4;
          return;
        }
        if((&p[ii]) == (&q[ii]))
        {
          a_didTestPass = false;
          a_errorCode = 5;
          return;
        }
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
        if(p0[ii] != 0)
        {
          a_didTestPass = false;
          a_errorCode = 6;
          return;
        }
        if(p1[ii] != 1)
        {
          a_didTestPass = false;
          a_errorCode = 7;
          return;
        }
        if(p2[ii] != 17)
        {
          a_didTestPass = false;
          a_errorCode = 8;
          return;
        }
        
        if ( (ii == 0) && (p3[ii] != 1))
        {
          a_didTestPass = false;
          a_errorCode = 9;
          return;
        }
        if ( (ii != 0) && (p3[ii] != 0))
        {
          a_didTestPass = false;
          a_errorCode = 10;
          return;
        }
        if ( (ii == DIM-1) && (p4[ii] != 17))
        {
          a_didTestPass = false;
          a_errorCode = 11;
          return;
        }
        if ( (ii != DIM-1) && (p4[ii] != 0))
        {
          a_didTestPass = false;
          a_errorCode = 12;
          return;
        }
      }
    }
    //Accessor Methods
    {
      Point p(1,2,3,4,5,6);
      for (int ii = 0; ii < DIM; ii++)
      {
        if(p[ii] != ii+1)
        {
          a_didTestPass = false;
          a_errorCode = 13;
          return;
        }
      }
      const Point q(1,2,3,4,5,6);
      for (int ii = 0; ii < DIM; ii++)
      {
        if(q[ii] != ii+1)
        {
          a_didTestPass = false;
          a_errorCode = 14;
          return;
        }

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
        if(p2[ii] != (p0[ii] + p1[ii]))
        {
          a_didTestPass = false;
          a_errorCode = 15;
          return;
        }
      }
    
      p2 = p1 - p0;
      for (int ii = 0; ii < DIM; ii++)
      {
        if(p2[ii] != (p1[ii] - p0[ii]))
        {
          a_didTestPass = false;
          a_errorCode = 16;
          return;
        }
      }

      p2 = p1*p0;
      for (int ii = 0; ii < DIM; ii++)
      {
        if(p2[ii] != (p1[ii]*p0[ii]))
        {
          a_didTestPass = false;
          a_errorCode = 17;
          return;
        }
      }

      p2 = p1/p0;
      for (int ii = 0; ii < DIM; ii++)
      {
        if(p2[ii] != (p1[ii]/p0[ii]))
        {
          a_didTestPass = false;
          a_errorCode = 18;
          return;
        }
      }

      p2 = p1/17;
      for (int ii = 0; ii < DIM; ii++)
      {
        if(p2[ii] != (p1[ii]/17))
        {
          a_didTestPass = false;
          a_errorCode = 18;
          return;
        }
      }

      p1 = Point::Ones(2); 
      p2 = p0 % p1;
      for (int ii = 0; ii < DIM; ii++)
      {
        if(p2[ii] != (p0[ii] % p1[ii]))
        {
          a_didTestPass = false;
          a_errorCode = 19;
          return;
        }
      }
    }  
  }
}
