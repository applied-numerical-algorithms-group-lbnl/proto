blt_add_library(NAME common
                HEADERS
                AMRRK4.H
                AMRSolver_FASMultigrid.H
                LevelSolver_FASMultigrid.H
                BoxOp_Laplace.H
                BoxOp_Template.H
                BoxOp_Euler.H
                LevelRK4.H
                RK4.H
                InputParser.H
                CommonTemplates.H
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/implem
                DEPENDS_ON Headers_AMR)
