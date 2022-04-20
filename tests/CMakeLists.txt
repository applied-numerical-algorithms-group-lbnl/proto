blt_add_executable(NAME BoxTests SOURCES BoxTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest)
blt_add_test(NAME BoxTests COMMAND BoxTests)
blt_add_executable(NAME BoxDataTests SOURCES BoxDataTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest)
blt_add_test(NAME BoxDataTests COMMAND BoxDataTests)
