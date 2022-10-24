blt_add_executable(NAME LevelEuler SOURCES main.cpp
    DEPENDS_ON Headers_Base common ${LIB_DEP})
if(ENABLE_HDF5 AND DEFINED ENV{HDF5_DIR})
    target_link_directories(LevelEuler PRIVATE $ENV{HDF5_DIR}/lib/)
    target_link_libraries(LevelEuler PRIVATE hdf5)
endif()
