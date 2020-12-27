
if ((Python_Development_FOUND) AND (USE_BRICK))
    include_directories(${Python_INCLUDE_DIRS})
    add_compile_definitions(BRICK_CODEGEN="$ENV{VecScatter_MODULE}" BRICK_COMPILER="${CMAKE_CXX_COMPILER}"
            BRICK_FLAGS="${CMAKE_CXX_FLAGS} -O2")

    set(SRCS jit.cpp)
    FORCUDA(${SRCS})

    add_library(brickJIT SHARED ${SRCS})

    set(SRCS main.cpp)
    FORCUDA(${SRCS})
    add_executable(brickJITexample ${SRCS})

    target_link_libraries(brickJITexample brick brickhelper ${EXTRALIBS})
endif()
