#include "Proto.H"
#include "Proto_Timer.H"
#include "UnitTestFunctions.H"
#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

#define _USE_MATH_DEFINES

using namespace Proto;
using namespace std;


int main(int argc, char** argv)
{
  //have to do this to get a time table
  PR_TIMER_SETFILE("proto.time.table");

  bool passed = false;
  int  errorCode = 0;

  prototest::pointTest(errorCode, passed);
  prototest::printTestMessage(string("Class Point   "), errorCode, passed);

  prototest::bxTest(errorCode, passed);
  prototest::printTestMessage(string("Class Bx      "), errorCode, passed);

  prototest::boxdataTest(errorCode, passed);
  prototest::printTestMessage(string("Class BoxData "), errorCode, passed);

  prototest::stencilTest(errorCode, passed);
  prototest::printTestMessage(string("Class Stencil "), errorCode, passed);

  prototest::interpTest(errorCode, passed);
  prototest::printTestMessage(string("Interp Stencil "), errorCode, passed);

  prototest::memtypeTest(errorCode, passed);
  prototest::printTestMessage(string("MemType "), errorCode, passed);


  PR_TIMER_REPORT();

} 

