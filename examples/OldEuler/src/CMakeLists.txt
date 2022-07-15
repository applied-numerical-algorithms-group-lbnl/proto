blt_add_library(NAME OldEuler_Src HEADERS EulerLevelDataRK4.H EulerOp.H
    SOURCES EulerLevelDataRK4.cpp EulerOp.cpp DEPENDS_ON Headers_Base common ${LIB_DEP}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
