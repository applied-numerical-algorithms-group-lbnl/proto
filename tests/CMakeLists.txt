blt_add_executable(NAME BoxTests SOURCES BoxTests.cpp
    DEPENDS_ON UnitTestFunctions ${LIB_DEP} gtest)
blt_add_test(NAME BoxTests COMMAND BoxTests)
