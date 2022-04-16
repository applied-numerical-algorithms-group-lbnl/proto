blt_add_library(NAME Headers_Base
                HEADERS  Proto.H
                Proto_Box.H
                Proto_BoxData.H
                Proto_DebugHooks.H
                Proto_MemType.H
                Proto_PAssert.H
                Proto_Point.H
                Proto_PowerItoI.H
                Proto_RK4.H
                Proto_Reduction.H
                Proto_Stencil.H
                Proto_Timer.H
                Proto_VisitWriter.H
                Proto_WriteBoxData.H
                Proto_gpu.H
                Proto_cuda.H
                Proto_HDF5.H
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/implem
                DEPENDS_ON proto_compiler_flags)
