# set(SRCS EulerOp.cpp GodunovAdvectionOp.cpp SGMultigrid.cpp)
set(SRCS GodunovAdvectionOp.cpp SGMultigrid.cpp)

FORCUDA(${SRCS})
add_library(common ${SRCS})
