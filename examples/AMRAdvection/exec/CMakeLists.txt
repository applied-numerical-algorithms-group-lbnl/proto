blt_add_executable(NAME AMRAdvection SOURCES main.cpp
    DEPENDS_ON AMRAdvection_Src Headers_AMR common ${LIB_DEP})
if(ENABLE_HDF5 AND DEFINED ENV{HDF5_DIR})
    target_link_directories(AMRAdvection PRIVATE $ENV{HDF5_DIR}/lib/)
    target_link_libraries(AMRAdvection PRIVATE hdf5)
endif()
