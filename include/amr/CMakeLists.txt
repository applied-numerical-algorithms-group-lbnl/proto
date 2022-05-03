blt_add_library(NAME Headers_AMR
                HEADERS ProtoAMR.H
                Proto_FinitePointSet.H
                Proto_LevelFluxRegister.H
                Proto_AMRGrid.H
                Proto_AMRData.H
                Proto_BoxOp.H
                Proto_LevelOp.H
                Proto_AMROp.H
                Proto_AMRSolver.H
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/implem ${CMAKE_CURRENT_SOURCE_DIR}/..
                DEPENDS_ON proto_compiler_flags Headers_Base)
