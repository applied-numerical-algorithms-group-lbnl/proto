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
if(ENABLE_CUDA)
    target_link_libraries(Headers_Base INTERFACE nvToolsExt)
endif()
if(ENABLE_HIP)
    target_include_directories(Headers_Base INTERFACE $ENV{ROCM_PATH}/roctracer/include/)
    target_link_directories(Headers_Base INTERFACE $ENV{ROCM_PATH}/lib)
    target_link_libraries(Headers_Base INTERFACE roctx64)
endif()
if(ENABLE_OPENMP)
    target_link_libraries(Headers_Base INTERFACE OpenMP::OpenMP_CXX)
endif()
if(ENABLE_HDF5)
    target_link_libraries(Headers_Base INTERFACE HDF5::HDF5)
endif()
