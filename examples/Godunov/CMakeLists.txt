set(SRCS exec/Godunov.cpp)
FORCUDA(${SRCS})

add_executable(godunov ${SRCS})
target_link_libraries(godunov common)
