#ifndef _PROTO_UNIT_TESTFunctions_
#define _PROTO_UNIT_TESTFunctions_

#include "Proto.H"

using std::string;
///
namespace prototest
{
  ///function to pring pass or fail and error code if failed
  void printTestMessage(const string& a_testName, const int & a_errorCode, bool a_didTestPass);

  ///function to test Proto::Point
  void pointTest(int& a_errorCode, bool & a_didTestPass);
}

#endif //end include guard
