blt_add_library(NAME UnitTestFunctions HEADERS UnitTestFunctions.H
                SOURCES UnitTestFunctions.cpp
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}
                DEPENDS_ON Headers_Base ${LIB_DEP})
