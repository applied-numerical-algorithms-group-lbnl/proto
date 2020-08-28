if (USE_BRICK)
    VSTARGET(MULTIGRID src/LevelMultigrid.cpp src/LevelMultigrid-out.cpp)

    set(MG_SRC ${VSTARGET_MULTIGRID_OUTPUT})
else()
    set(MG_SRC src/LevelMultigrid.cpp)
endif()

set(SRCS exec/multigridTest.cpp ${MG_SRC})
message(STATUS ${SRCS})
FORCUDA(${SRCS})

include_directories(src)
add_executable(multigridleveldata ${SRCS})

if (USE_BRICK)
    target_link_libraries(multigridleveldata brick brickhelper)
endif()
