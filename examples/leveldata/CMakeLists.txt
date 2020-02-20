set(SRCS exec/streams.cpp)
FORCUDA(${SRCS})

add_executable(leveldata ${SRCS})
target_link_libraries(leveldata common)
