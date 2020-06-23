cmake_minimum_required(VERSION 3.13)
project(proto
        LANGUAGES CXX)

cmake_policy(SET CMP0077 NEW)

include(CTest)
add_custom_target(BuildTests)
add_custom_target(check-all COMMAND ${CMAKE_CTEST_COMMAND} --progress)
add_dependencies(check-all BuildTests)

# Usage:
#  - mkdir build && cd build
#  - cmake ..
# -DXX=YY to specify options

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(IgnoreMe "${CMAKE_C_COMPILER}")

include(CheckLanguage)
include(CheckCXXCompilerFlag)

# define target Coverage for GNU compilers
set(CMAKE_CXX_FLAGS_COVERAGE "-g -O0 -fprofile-arcs -ftest-coverage" CACHE STRING
        "Flags used by the C++ compiler during coverage builds."
        FORCE)
set(CMAKE_C_FLAGS_COVERAGE "-g -O0 -fprofile-arcs -ftest-coverage" CACHE STRING
        "Flags used by the C compiler during coverage builds."
        FORCE)
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE
        "-fprofile-arcs -ftest-coverage" CACHE STRING
        "Flags used for linking binaries during coverage builds."
        FORCE)

mark_as_advanced(
        CMAKE_CXX_FLAGS_COVERAGE
        CMAKE_C_FLAGS_COVERAGE
        CMAKE_EXE_LINKER_FLAGS_COVERAGE
        CMAKE_SHARED_LINKER_FLAGS_COVERAGE)
# Update the documentation string of CMAKE_BUILD_TYPE for GUIs
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel Coverage."
        FORCE)

# Define coverage as a build target to generate HTML report
if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
    find_program(GCOV_PATH gcov)
    find_program(GCOVR_PATH gcovr)
    if((NOT GCOV_PATH) OR (NOT GCOVR_PATH))
        message(WARNING "gcov(r) not found! Unable to generate coverage report")
    else()
        add_custom_target(coverage
                # Create folder
                COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/coverage

                # Running gcovr
                COMMAND ${GCOVR_PATH} --html --html-details
                -r ${CMAKE_SOURCE_DIR}
                -f ${CMAKE_SOURCE_DIR}/include
                --object-directory=${PROJECT_BINARY_DIR}
                -o coverage/index.html

                BYPRODUCTS ${PROJECT_BINARY_DIR}/coverage  # report directory
                WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                DEPENDS check-all
                VERBATIM # Protect arguments to commands
                COMMENT "Running gcovr to produce HTML code coverage report."
        )

        add_custom_command(TARGET coverage POST_BUILD
                COMMAND ;
                COMMENT "Open ./coverage/index.html in your browser to view the coverage report."
        )
    endif()
endif()

# Use OpenMP
find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DTHRUST_DEBUG")

check_cxx_compiler_flag(-march=native HAS_MARCH)
if (HAS_MARCH)
    message(STATUS "Using march=native")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
else ()
    # Work around Summit only have mcpu flag
    check_cxx_compiler_flag(-mcpu=native HAS_MCPU)
    if (HAS_MCPU)
        message(STATUS "Using mcpu=native")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=native")
    endif ()
endif ()

# Option - DIM (2/3)
set(DIM 3 CACHE STRING "Number of dimensions")
set_property(CACHE DIM PROPERTY STRINGS 2 3)
if (DIM)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDIM=${DIM}")
endif ()

# Option - MEM_CHECK (ON/OFF)
option(MEM_CHECK "Compile with memcheck" OFF)
if (MEM_CHECK)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPROTO_MEM_CHECK")
endif ()

# Option - TIMER (ON/OFF)
option(TIMER "Compiler with timing" ON)
if (NOT TIMER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPR_TURN_OFF_TIMERS")
endif ()

# Option - USE_CUDA (ON/OFF)
option(USE_CUDA "This build uses CUDA" OFF)
if (USE_CUDA)
    check_language(CUDA)

    if (CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})
        set(CMAKE_CUDA_SEPARABLE_COMPILATION OFF)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPROTO_CUDA=1")
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler \"${CMAKE_CXX_FLAGS}\"")
        # This is a hack for CMake to compile .cpp files using nvcc
        function(FORCUDA)
            set_source_files_properties(${ARGV} PROPERTIES LANGUAGE CUDA)
        endfunction()
    endif ()
endif ()

if (NOT (COMMAND FORCUDA))
    function(FORCUDA)
    endfunction()
endif ()

# Option - USE_BRICK (ON/OFF)
option(USE_BRICK "Enable building with brick" ON)
if (USE_BRICK)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPROTO_BRICK")
    set(BRICK_PROJECT_PATH ${PROJECT_SOURCE_DIR}/external/bricklib)
    include(${BRICK_PROJECT_PATH}/cmake/VecScatter.cmake)
    set(BRICK_BUILD_TEST OFF)
    add_subdirectory(
            ${PROJECT_SOURCE_DIR}/external/bricklib
    )
    include_directories(${BRICK_PROJECT_PATH}/include)
endif ()

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/EBProto/include
)

add_subdirectory(examples)
