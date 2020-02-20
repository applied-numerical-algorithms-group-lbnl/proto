set(SRCS exec/applyKernel.cpp)
FORCUDA(${SRCS})

add_executable(applykernel ${SRCS})
target_link_libraries(applykernel common)
