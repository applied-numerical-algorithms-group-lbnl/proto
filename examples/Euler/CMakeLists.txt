set(SRCS exec/Euler.cpp src/EulerRK4.cpp)
FORCUDA(${SRCS})

include_directories(src)
add_executable(euler ${SRCS})
target_link_libraries(euler common)
