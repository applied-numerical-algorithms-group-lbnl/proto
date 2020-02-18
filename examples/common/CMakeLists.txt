set(SRCS EulerOp.cpp GodunovAdvectionOp.cpp SGMultigrid.cpp)

FORCUDA(${SRCS})
add_library(common ${SRCS})
