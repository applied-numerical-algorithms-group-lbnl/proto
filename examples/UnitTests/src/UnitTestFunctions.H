#ifndef _PROTO_UNIT_TESTFunctions_
#define _PROTO_UNIT_TESTFunctions_

#include "Proto.H"

using std::string;
///
namespace prototest
{
  ///for debuggging
  void breakHere(int a_errorCode);

  /// returns true if test passed, if test fails, sets error code to test number
  bool UNIT_TEST(bool a_test, int& a_errorCode, int testnumber);

  ///function to pring pass or fail and error code if failed
  void printTestMessage(const string& a_testName, const int & a_errorCode, bool a_didTestPass);

  ///function to test Proto::Point
  void pointTest(int& a_errorCode, bool & a_didTestPass);


  ///function to test Proto::Box
  void bxTest(int& a_errorCode, bool & a_didTestPass);


  ///function to test Proto::BoxData
  void boxdataTest(int& a_errorCode, bool & a_didTestPass);


  ///function to test Proto::BoxData
  void stencilTest(int& a_errorCode, bool & a_didTestPass);

  ///function to test Proto::BoxData
  void interpTest(int& a_errorCode, bool & a_didTestPass);

   ///function to test Proto::MEMTYPE
  void memtypeTest(int& a_errorCode, bool & a_didTestPass);
}

#endif //end include guard
