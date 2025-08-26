list(APPEND COMMON_HEADERS
    InputParser.H
    BoxOp_Laplace.H
    BoxOp_Template.H
    BoxOp_Euler.H)
blt_add_library(NAME common
  HEADERS ${COMMON_HEADERS}
  INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} 
  DEPENDS_ON Headers_AMR ${LIB_DEP})
