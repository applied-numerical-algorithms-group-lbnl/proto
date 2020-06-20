add_unittest(UT_BoxData cInterval.cpp boxData.cpp)
if (USE_BRICK)
    target_link_libraries(UT_BoxData brick brickhelper)
endif()
