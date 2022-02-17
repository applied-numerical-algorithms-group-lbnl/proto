blt_add_library(NAME LevelMultigrid HEADERS LevelMultigrid.H SOURCES LevelMultigrid.cpp
    DEPENDS_ON Headers ${LIB_DEP} INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
