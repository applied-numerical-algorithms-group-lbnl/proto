if(ENABLE_HDF5)
    blt_add_executable(NAME HDF5 SOURCES HDF5.cpp DEPENDS_ON Headers ${LIB_DEP})
    target_link_libraries(HDF5 PRIVATE HDF5::HDF5)
endif()
