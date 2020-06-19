include_directories(${GTEST_INCLUDE_DIRS})

function(add_unittest testname)
    add_executable(${testname} EXCLUDE_FROM_ALL ${ARGN})
    gtest_add_tests(TARGET ${testname})
    target_link_libraries(${testname} ${UTEST_MODULES} ${GTEST_BOTH_LIBRARIES})
    add_dependencies(BuildTests ${testname})
endfunction()

add_subdirectory(basic)
add_subdirectory(Point)
add_subdirectory(Box)
