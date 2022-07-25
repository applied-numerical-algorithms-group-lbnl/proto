blt_add_library(NAME Particle HEADERS particle.hpp
    DEPENDS_ON Headers common ${LIB_DEP}
    INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
