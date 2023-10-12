blt_add_executable(NAME FourthOrder1DWENO SOURCES FourthOrder1DWENO.cpp
    DEPENDS_ON AdvectionRK4 ${LIB_DEP})
