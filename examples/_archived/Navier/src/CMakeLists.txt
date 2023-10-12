blt_add_library(NAME BCG_Integrator HEADERS BCG_Integrator.H SOURCES BCG_Integrator.cpp
    DEPENDS_ON common Headers ${LIB_DEP} INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
