blt_add_executable(NAME EulerLevelData SOURCES EulerTest.cpp
    DEPENDS_ON EulerLevelDataRK4 Headers ${LIB_DEP})
