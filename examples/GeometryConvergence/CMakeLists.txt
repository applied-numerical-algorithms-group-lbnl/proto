set(SRCS exec/GeometryConvergence.cpp)
FORCUDA(${SRCS})

add_executable(geometryConvergence ${SRCS})
target_link_libraries(geometryConvergence common)
