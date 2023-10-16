list(APPEND COMMON_HEADERS
    LevelSolver_FASMultigrid.H
    RK4.H
    InputParser.H
    CommonTemplates.H 
    LevelRK4.H
    BoxOp_Laplace.H
    BoxOp_Template.H
    BoxOp_Euler.H)
if(AMR)
    list(APPEND COMMON_HEADERS
        AMRRK4.H
        AMRSolver_FASMultigrid.H)
endif()
blt_add_library(NAME common
  HEADERS ${COMMON_HEADERS}
  INCLUDES ${CMAKE_CURRENT_SOURCE_DIR} 
  DEPENDS_ON Headers_AMR ${LIB_DEP})
