blt_add_library(NAME SolidRot HEADERS SRK4.H particle.hpp
    DEPENDS_ON Headers common ${LIB_DEP}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
