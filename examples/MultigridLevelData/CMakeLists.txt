set(SRCS exec/multigridTest.cpp src/LevelMultigrid.cpp)
FORCUDA(${SRCS})

include_directories(src)
add_executable(multigridleveldata ${SRCS})

if (USE_BRICK)
    target_link_libraries(multigridleveldata brick brickhelper)
endif()
