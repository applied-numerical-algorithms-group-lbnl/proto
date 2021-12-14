blt_add_executable(NAME BoxDataUnit SOURCES BoxDataUnit.cpp
    DEPENDS_ON UnitTestFunctions ${LIB_DEP} gtest)
blt_add_test(NAME BoxDataUnit COMMAND BoxDataUnit)
