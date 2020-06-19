if (USE_BRICK)
    VSTARGET(N3MPI_SHMEM exec/streams.cpp exec/streams-out.cpp)

    set(SRCS ${VSTARGET_N3MPI_SHMEM_OUTPUT})
else()
    set(SRCS exec/streams.cpp)
endif()

FORCUDA(${SRCS})

add_executable(laplacian ${SRCS})

if (USE_BRICK)
    target_link_libraries(laplacian brick brickhelper)
endif()

add_test(NAME laplacian COMMAND laplacian)
add_dependencies(BuildTests laplacian)
