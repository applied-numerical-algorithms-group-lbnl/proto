blt_add_library(NAME common HEADERS CommonTemplates.H EulerOp.H GodunovAdvectionOp.H SGMultigrid.H
                SOURCES EulerOp.cpp GodunovAdvectionOp.cpp SGMultigrid.cpp
                DEPENDS_ON Headers ${LIB_DEP} INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
