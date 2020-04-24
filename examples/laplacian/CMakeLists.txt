VSTARGET(N3MPI_SHMEM exec/streams.cpp exec/streams-out.cpp)

set(SRCS ${VSTARGET_N3MPI_SHMEM_OUTPUT})
FORCUDA(${SRCS})

add_executable(laplacian ${SRCS})
target_link_libraries(laplacian common)
