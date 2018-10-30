#include "Proto.H"
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
  bool passed = false;
  int  errorCode = 0;

  prototest::pointTest(errorCode, passed);
  prototest::printTestMessage(string("Class Point   "), errorCode, passed);

  prototest::bxTest(errorCode, passed);
  prototest::printTestMessage(string("Class Bx      "), errorCode, passed);

  prototest::boxdataTest(errorCode, passed);
  prototest::printTestMessage(string("Class BoxData "), errorCode, passed);

} // end main

