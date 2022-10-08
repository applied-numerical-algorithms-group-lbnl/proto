blt_add_executable(NAME FASMultigrid SOURCES main.cpp
    DEPENDS_ON Headers_AMR common ${LIB_DEP})
if(ENABLE_HDF5 AND DEFINED ENV{HDF5_DIR})
    target_link_directories(FASMultigrid PRIVATE $ENV{HDF5_DIR}/lib/)
    target_link_libraries(FASMultigrid PRIVATE hdf5)
endif()
