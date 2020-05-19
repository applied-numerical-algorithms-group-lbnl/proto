if (USE_BRICK)
    VSTARGET(LAPLACIAN exec/streams.cpp exec/streams-out.cpp)

    set(SRCS ${VSTARGET_LAPLACIAN_OUTPUT})
else()
    set(SRCS exec/streams.cpp)
endif()

FORCUDA(${SRCS})

add_executable(laplacian ${SRCS})

if (USE_BRICK)
    target_link_libraries(laplacian common brick brickhelper)
else()
    target_link_libraries(laplacian common)
endif()
