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
if(ENABLE_CUDA)
    target_link_libraries(Headers_AMR INTERFACE nvToolsExt)
endif()
if(ENABLE_HIP)
    target_include_directories(Headers_AMR INTERFACE $ENV{ROCM_PATH}/roctracer/include/)
    target_link_directories(Headers_AMR INTERFACE $ENV{ROCM_PATH}/lib)
    target_link_libraries(Headers_AMR INTERFACE roctx64)
endif()
if(ENABLE_OPENMP)
    target_link_libraries(Headers_AMR INTERFACE OpenMP::OpenMP_CXX)
endif()
if(ENABLE_HDF5)
    target_link_libraries(Headers_AMR INTERFACE HDF5::HDF5)
endif()
