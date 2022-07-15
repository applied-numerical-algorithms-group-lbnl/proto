blt_add_library(NAME common
                HEADERS
                AMRSolver_FASMultigrid.H
                LevelSolver_FASMultigrid.H
                BoxOp_Laplace.H
                InputParser.H
                CommonTemplates.H
                BoxOp_Euler.H
                LevelRK4.H
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/implem
                DEPENDS_ON proto_compiler_flags Headers_AMR)
