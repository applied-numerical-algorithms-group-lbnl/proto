blt_add_library(NAME AMRAdvection_Src HEADERS BoxOp_Advection.H
    DEPENDS_ON Headers_AMR common ${LIB_DEP} INCLUDES ${CMAKE_CURRENT_SOURCE_DIR})
