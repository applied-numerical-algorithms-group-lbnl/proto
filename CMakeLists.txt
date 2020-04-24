cmake_minimum_required(VERSION 3.13)
project(proto
        LANGUAGES CXX)

cmake_policy(SET CMP0077 NEW)

# Usage:
#  - mkdir build && cd build
#  - cmake ..
# -DXX=YY to specify options

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CheckLanguage)

# Use OpenMP
#find_package(OpenMP REQUIRED)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DTHRUST_DEBUG")

# Option - DIM (2/3)
set(DIM 2 CACHE STRING "Number of dimensions")
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

set(BRICK_PROJECT_PATH ${PROJECT_SOURCE_DIR}/external/bricklib)
include(${BRICK_PROJECT_PATH}/cmake/VecScatter.cmake)
set(BRICK_BUILD_TEST OFF)
add_subdirectory(
        ${PROJECT_SOURCE_DIR}/external/bricklib
)

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/EBProto/include
        ${BRICK_PROJECT_PATH}/include
)

add_subdirectory(examples)
