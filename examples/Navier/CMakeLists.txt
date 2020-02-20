set(SRCS exec/Navier.cpp src/BCG_Integrator.cpp)
FORCUDA(${SRCS})

include_directories(src)
add_executable(navier ${SRCS})
target_link_libraries(navier common)
