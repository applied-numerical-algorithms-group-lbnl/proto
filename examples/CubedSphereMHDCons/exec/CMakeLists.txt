blt_add_executable(NAME CubedSphereMHDCons SOURCES main.cpp
    DEPENDS_ON CubedSphereMHDCons_Src ${LIB_DEP})
target_compile_definitions(CubedSphereMHDCons PRIVATE
        CUBED_SPHERE_SHELL_R0=0.1 CUBED_SPHERE_SHELL_R1=2.0
        CUBED_SPHERE_SHELL_RADIAL_COORD=0 CUBED_SPHERE_RADIAL_STRETCH=1.5)
