#ifndef _PROTO_UNIT_TEST_
#define _PROTO_UNIT_TEST_

static int VERBO;
static int s_fillWidth = 50;
static bool PASS = true;
static bool ALLPASS = true;
static int TEST;

/// Macro Based Unit Testing API
/** Basic Structure:

    #include "UnitTest.H"

    @code
    int main(int argc, char** argv)
    {
      // Set verbosity and test num
      if (argc == 3)
      {
        TEST = atoi(argv[1]);
        VERBO = atoi(argv[2]);
      }
      else
      {
        TEST = 0;
        VERBO = 0;
      }

      // Begin a test suite
      // If the first argument agrees with TEST, this suite will be run
      // The string argument will be printed
      BEGIN_TEST_SUITE(0, "Test Suite Zero")
      BEGIN_TEST("First Test")
      // code here
      
      // print output contingent on verbosity. Code between __OUT and OUT__
      // will run if verbosity is equal to or high than the argument of __OUT
      __OUT(2) std::cout << "Verbosity is at least 2!" << std::endl; OUT__

      // Between BEGIN_TEST and END_TEST, call UNIT_TEST any number of times.
      // It is useful to place a full statement inside UNIT_TEST because
      // that will be the text that is printed when it fails.
      UNIT_TEST(<statement that returns a bool>)
      UNIT_TEST(<another statement>) //etc
      
      // Done Testing
      END_TEST

      BEGIN_TEST("Second Test")
      // ...
      END_TEST

      // ...

      END_TEST_SUITE
      @endcode
    }

    verbosity defaults:
      VERBO = 0:  Test Suite level events. Only prints the name of the test suite and whether all tests were successful
      VERBO = 1:  Individual test level events. Prints the name of each test and whether it was successful
      VERBO > 1:  Unused by default. Recommended for sub-test level events such as printing intermediate output

*/

/// Dummy function for placing breakpoints.
int _breakPoint()
{
  int bpt = 1; //put a break point here
  return bpt + 1;
}

#define BEGIN_TEST_SUITE(test_num, test_name) \
    if (TEST == test_num) {                                               \
    cout << endl << "Test " << TEST << ": " << test_name << endl;         \
    if (VERBO > 0)                                                  \
    {                                                                     \
        cout << setfill('=') << setw(s_fillWidth) << " " <<  endl;        \
    }                                                                     \
    PASS = true;                                                          \
    ALLPASS = true;                                                       \
    {

#define BEGIN_TEST(test_name)                                             \
    if (VERBO > 0)                                                  \
    {                                                                     \
        cout << endl << "Testing " << test_name << endl;                  \
        cout << setfill('-') << setw(s_fillWidth) << " " <<  endl;        \
    }                                                                     \
    PASS = true;                                                          \
    {

#define END_TEST()                                                        \
    }                                                                     \
    ALLPASS &= PASS;                                                      \
    if (VERBO > 0)                                                  \
    {                                                                     \
        if (PASS) {cout << endl << "TEST PASSED" << endl;}                \
        else {cout << endl << "TEST FAILED" << endl;}                     \
    }                                                                     \
    
#define END_TEST_SUITE()                                                  \
    }                                                                     \
    if (VERBO > 0)                                                  \
    {                                                                     \
        cout << setfill('-') << setw(s_fillWidth) << " " <<  endl;        \
    }                                                                     \
    if (ALLPASS) {cout << "All Tests Successful" << endl;}                \
    else {cout << "Some Tests Failed" << endl;}                           \
    if (VERBO > 0)                                                  \
    {                                                                     \
        cout << setfill('=') << setw(s_fillWidth) << " " <<  endl;        \
    }                                                                     \
    }                                                                     \

#define OMIT_TEST(test_name, msg)                                         \
    cout << endl << setfill('-') << setw(s_fillWidth) << " " <<  endl;    \
    cout << "TEST OMITTED: " << test_name << endl;                        \
    cout << "\t" << msg << endl;                                          \

#define UNIT_TEST(x)                                                      \
    PASS &= x;                                                            \
    if (!(x)) \
    {\
      cout << "Line " << __LINE__ << " -- Unit Test Failed: " << #x << endl;\
      _breakPoint();\
    }

#define __OUT(level) \
    if (VERBO >= int(level)) {\


#define OUT__ }
#endif //end include guard
