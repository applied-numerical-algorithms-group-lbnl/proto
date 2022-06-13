blt_add_library(NAME common
                HEADERS
                AMRRK4.H
                AMRSolver_FASMultigrid.H
                LevelSolver_FASMultigrid.H
                BoxOp_Laplace.H
                BoxOp_Template.H
                InputParser.H
                TestFunc.H
                amr_subcycle_v1/AMRSubcycleExplicit.H
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} 
                DEPENDS_ON Headers_AMR)
