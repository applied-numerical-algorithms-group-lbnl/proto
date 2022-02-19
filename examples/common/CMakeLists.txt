blt_add_library(NAME common
                HEADERS
                AMRSolver_FASMultigrid.H
                LevelSolver_FASMultigrid.H
                BoxOp_Laplace.H
                INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/implem
                DEPENDS_ON proto_compiler_flags Headers_AMR)
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
