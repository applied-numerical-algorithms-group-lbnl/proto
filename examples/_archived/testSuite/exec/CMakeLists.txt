blt_add_executable(NAME testSuite SOURCES testSuite.cpp
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS_ON Headers ${LIB_DEP})
