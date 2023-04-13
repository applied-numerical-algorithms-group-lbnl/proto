blt_add_executable(NAME AMRFAS SOURCES main.cpp
    DEPENDS_ON Headers_AMR common ${LIB_DEP})
if(DEFINED ENV{OLCF_OPENBLAS_ROOT})
    target_link_directories(AMRFAS PRIVATE $ENV{OLCF_OPENBLAS_ROOT}/lib/)
    target_link_libraries(AMRFAS PRIVATE openblas)
endif()
