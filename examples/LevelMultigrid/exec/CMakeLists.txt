blt_add_executable(NAME LevelMultigrid SOURCES main.cpp
    DEPENDS_ON LevelMultigrid_Src ${LIB_DEP})
if(ENABLE_HDF5 AND DEFINED ENV{HDF5_DIR})
    target_link_directories(LevelMultigrid PRIVATE $ENV{HDF5_DIR}/lib/)
    target_link_libraries(LevelMultigrid PRIVATE hdf5)
endif()
