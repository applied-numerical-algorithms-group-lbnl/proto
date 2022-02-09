#configure_file(Proto_Config.H.in Proto_Config.H)
blt_add_library(NAME Headers
                HEADERS  Proto.H
                Proto_Box.H
                Proto_BoxData.H
                Proto_Config.H.in
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
    target_link_libraries(Headers INTERFACE nvToolsExt)
endif()
if(ENABLE_HIP)
    target_include_directories(Headers INTERFACE $ENV{ROCM_PATH}/roctracer/include/)
    target_link_directories(Headers INTERFACE $ENV{ROCM_PATH}/lib)
    target_link_libraries(Headers INTERFACE roctx64)
endif()
if(ENABLE_OPENMP)
    target_link_libraries(Headers INTERFACE OpenMP::OpenMP_CXX)
endif()
if(ENABLE_HDF5)
    target_link_libraries(Headers INTERFACE HDF5::HDF5)
endif()
