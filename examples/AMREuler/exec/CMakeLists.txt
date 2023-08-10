blt_add_executable(NAME AMREuler SOURCES main.cpp
    DEPENDS_ON Headers_AMR common ${LIB_DEP})
if(ENABLE_HDF5 AND HDF5_FOUND)
	target_link_libraries(AMREuler PUBLIC HDF5::HDF5)
endif()
