set(SRCS exec/PoissonStencil.cpp)
FORCUDA(${SRCS})

add_executable(poisson ${SRCS})
target_link_libraries(poisson common)
