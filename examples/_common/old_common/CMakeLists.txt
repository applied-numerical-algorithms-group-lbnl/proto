blt_add_library(NAME common HEADERS CommonTemplates.H EulerOp.H GodunovAdvectionOp.H SGMultigrid.H
                AMRSolver_FASMultigrid.H LevelSolver_FASMultigrid.H BoxOp_Laplace.H
                SOURCES EulerOp.cpp GodunovAdvectionOp.cpp SGMultigrid.cpp
                DEPENDS_ON Headers_AMR ${LIB_DEP} INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
