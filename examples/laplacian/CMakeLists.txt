set(SRCS exec/streams.cpp)
FORCUDA(${SRCS})

add_executable(laplacian ${SRCS})
target_link_libraries(laplacian common)
