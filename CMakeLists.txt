message(STATUS "CMake started")

cmake_minimum_required(VERSION 3.2)

project(OPTICS-Clustering)


set(DEBUG_FLAGS "-g -O1 -pthread -pedantic")
set(RELEASE_FLAGS "-O3 -pthread -pedantic")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS ${RELEASE_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})
set(CMAKE_CONFIGURATION_TYPES Debug Release)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()

add_subdirectory(test)