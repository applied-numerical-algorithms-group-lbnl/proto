blt_add_executable(NAME AMRFAS SOURCES main.cpp 
    DEPENDS_ON Headers_AMR common ${LIB_DEP} LAPACK::LAPACK)
