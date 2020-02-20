set(SRCS exec/Multigrid.cpp)
FORCUDA(${SRCS})

add_executable(multigrid ${SRCS})
target_link_libraries(multigrid common)
