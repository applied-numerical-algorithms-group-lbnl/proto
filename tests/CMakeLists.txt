blt_add_executable(NAME BoxTests SOURCES BoxTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME BoxTests COMMAND BoxTests)
blt_add_executable(NAME ReductionTests SOURCES ReductionTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME ReductionTests COMMAND ReductionTests)
blt_add_executable(NAME BoxDataTests SOURCES BoxDataTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME BoxDataTests COMMAND BoxDataTests)
blt_add_executable(NAME forAllTests SOURCES forAllTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME forAllTests COMMAND forAllTests)
blt_add_executable(NAME StencilTests SOURCES StencilTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME StencilTests COMMAND StencilTests)
blt_add_executable(NAME InterpTests SOURCES InterpTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME InterpTests COMMAND InterpTests)
