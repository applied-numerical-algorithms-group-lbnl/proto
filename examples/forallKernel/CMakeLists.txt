set(SRCS exec/forallKernel.cpp)
FORCUDA(${SRCS})

add_executable(forallKernel ${SRCS})
target_link_libraries(forallKernel common)
