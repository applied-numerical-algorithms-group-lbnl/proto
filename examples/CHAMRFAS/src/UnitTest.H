#pragma once

static int _VERBO_;
static int s_fillWidth = 50;
static bool _PASS_ = true;
static bool _ALL_PASS_ = true;
static int _TEST_;

using std::cout;
using std::endl;

#define BEGIN_TEST_SUITE(test_num, test_name) \
    if (_TEST_ == test_num) {                                               \
    cout << endl << "Test " << _TEST_ << ": " << test_name << endl;         \
    if (_VERBO_ > 0)                                                  \
    {                                                                     \
        cout << setfill('=') << setw(s_fillWidth) << " " <<  endl;        \
    }                                                                     \
    _PASS_ = true;                                                          \
    _ALL_PASS_ = true;                                                       \
    {

#define BEGIN_TEST(test_name)                                             \
    if (_VERBO_ > 0)                                                  \
    {                                                                     \
        cout << endl << "Testing " << test_name << endl;                  \
        cout << setfill('-') << setw(s_fillWidth) << " " <<  endl;        \
    }                                                                     \
    _PASS_ = true;                                                          \
    {

#define END_TEST()                                                        \
    }                                                                     \
    _ALL_PASS_ &= _PASS_;                                                      \
    if (_VERBO_ > 0)                                                  \
    {                                                                     \
        if (_PASS_) {cout << endl << "TEST PASSED" << endl;}                \
        else {cout << endl << "TEST FAILED" << endl;}                     \
    }                                                                     \
    
#define END_TEST_SUITE()                                                  \
    }                                                                     \
    if (_VERBO_ > 0)                                                  \
    {                                                                     \
        cout << setfill('-') << setw(s_fillWidth) << " " <<  endl;        \
    }                                                                     \
    if (_ALL_PASS_) {cout << "All Tests Successful" << endl;}                \
    else {cout << "Some Tests Failed" << endl;}                           \
    if (_VERBO_ > 0)                                                  \
    {                                                                     \
        cout << setfill('=') << setw(s_fillWidth) << " " <<  endl;        \
    }                                                                     \
    }                                                                     \

#define OMIT_TEST(test_name, msg)                                         \
    cout << endl << setfill('-') << setw(s_fillWidth) << " " <<  endl;    \
    cout << "TEST OMITTED: " << test_name << endl;                        \
    cout << "\t" << msg << endl;                                          \

#define UNIT_TEST(x)                                                      \
    _PASS_ &= x;                                                            \
    if (_VERBO_ > 2){                                                       \
    if (!x) {cout << "Unit Test Failed: " << #x << endl;}\
    }\

void unitTest(bool success){_PASS_ &= success;}
