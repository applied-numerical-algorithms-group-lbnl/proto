blt_add_executable(NAME PointTests SOURCES PointTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME PointTests COMMAND PointTests)
blt_add_executable(NAME BoxTests SOURCES BoxTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME BoxTests COMMAND BoxTests)
blt_add_executable(NAME CoordPermutationTests SOURCES CoordPermutationTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME CoordPermutationTests COMMAND CoordPermutationTests)
blt_add_executable(NAME ReductionTests SOURCES ReductionTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME ReductionTests COMMAND ReductionTests)
blt_add_executable(NAME BoxDataTests SOURCES BoxDataTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME BoxDataTests COMMAND BoxDataTests)
blt_add_executable(NAME ForallTests SOURCES ForallTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME ForallTests COMMAND ForallTests)
blt_add_executable(NAME StencilTests SOURCES StencilTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME StencilTests COMMAND StencilTests)
blt_add_executable(NAME ArrayTests SOURCES ArrayTests.cpp
    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
blt_add_test(NAME ArrayTests COMMAND ArrayTests)
#blt_add_executable(NAME InterpTests SOURCES InterpTests.cpp
#    DEPENDS_ON Headers_Base ${LIB_DEP} gtest
#    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
#blt_add_test(NAME InterpTests COMMAND InterpTests)
if(OPS)
  blt_add_executable(NAME MatrixTests SOURCES MatrixTests.cpp
      DEPENDS_ON Headers_OPS ${LIB_DEP} gtest
      INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
  blt_add_test(NAME MatrixTests COMMAND MatrixTests)
endif()
if(ENABLE_ALL_TESTS)
    blt_add_executable(NAME DisjointBoxLayoutTests SOURCES DisjointBoxLayoutTests.cpp
        DEPENDS_ON Headers_Base ${LIB_DEP} gtest
        INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
    blt_add_test(NAME DisjointBoxLayoutTests COMMAND DisjointBoxLayoutTests)
    blt_add_executable(NAME LevelBoxDataTests SOURCES LevelBoxDataTests.cpp
        DEPENDS_ON Headers_Base ${LIB_DEP} gtest
        INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
    blt_add_test(NAME LevelBoxDataTests COMMAND LevelBoxDataTests)
    if (AMR)
        blt_add_executable(NAME AMRDataTests SOURCES AMRDataTests.cpp
            DEPENDS_ON Headers_AMR ${LIB_DEP} gtest
            INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
        blt_add_test(NAME AMRDataTests COMMAND AMRDataTests)
        blt_add_executable(NAME LevelFluxRegisterTests SOURCES LevelFluxRegisterTests.cpp
            DEPENDS_ON Headers_AMR ${LIB_DEP} gtest
            INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
        blt_add_test(NAME LevelFluxRegisterTests COMMAND LevelFluxRegisterTests)
    endif()
    if (MMB)
        blt_add_executable(NAME MBGraphTests SOURCES MBGraphTests.cpp
            DEPENDS_ON Headers_MMB ${LIB_DEP} gtest
            INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
        blt_add_test(NAME MBGraphTests COMMAND MBGraphTests)
        blt_add_executable(NAME MBProblemDomainTests SOURCES MBProblemDomainTests.cpp
            DEPENDS_ON Headers_MMB ${LIB_DEP} gtest
            INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
        blt_add_test(NAME MBProblemDomainTests COMMAND MBProblemDomainTests)
        blt_add_executable(NAME MBDisjointBoxLayoutTests SOURCES MBDisjointBoxLayoutTests.cpp
            DEPENDS_ON Headers_MMB ${LIB_DEP} gtest
            INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
        blt_add_test(NAME MBDisjointBoxLayoutTests COMMAND MBDisjointBoxLayoutTests)
        blt_add_executable(NAME MBLevelBoxDataTests SOURCES MBLevelBoxDataTests.cpp
            DEPENDS_ON Headers_MMB ${LIB_DEP} gtest
            INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
        blt_add_test(NAME MBLevelBoxDataTests COMMAND MBLevelBoxDataTests)
        blt_add_executable(NAME MBMapTests SOURCES MBMapTests.cpp
            DEPENDS_ON Headers_MMB ${LIB_DEP} gtest
            INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
        blt_add_test(NAME MBMapTests COMMAND MBMapTests)
    endif()
    if (OPS)
        blt_add_executable(NAME OperatorTests SOURCES OperatorTests.cpp
          DEPENDS_ON Headers_OPS ${LIB_DEP} gtest
          INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
        blt_add_test(NAME OperatorTests COMMAND OperatorTests)
    endif()
    if (MMB AND OPS)
        blt_add_executable(NAME MBInterpOpTests SOURCES MBInterpOpTests.cpp
            DEPENDS_ON Headers_MMB Headers_OPS ${LIB_DEP} gtest
            INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
        blt_add_test(NAME MBInterpOpTests COMMAND MBInterpOpTests)
    endif()
endif()
