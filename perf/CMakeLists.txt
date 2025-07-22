
blt_add_executable(NAME ReductionPerf SOURCES ReductionPerf.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME ReductionPerf COMMAND ReductionPerf)
blt_add_executable(NAME LaplaceStencilPerf SOURCES LaplaceStencilPerf.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME LaplaceStencilPerf COMMAND LaplaceStencilPerf)