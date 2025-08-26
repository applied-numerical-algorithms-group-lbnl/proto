blt_add_executable(NAME AMRPoisson SOURCES main.cpp
    DEPENDS_ON Headers_AMR common ${LIB_DEP})
if(DEFINED ENV{OLCF_OPENBLAS_ROOT})
    target_link_directories(AMRPoisson PRIVATE $ENV{OLCF_OPENBLAS_ROOT}/lib/)
    target_link_libraries(AMRPoisson PRIVATE openblas)
endif()
