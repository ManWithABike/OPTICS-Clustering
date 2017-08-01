message(STATUS "CMake started")

cmake_minimum_required(VERSION 2.8)

project(OPTICS-Clustering)

if(CMAKE_COMPILER_IS_GNUCXX)
#gcc needs additional linker-flags "-lX11" & "-pthread"
	set(DEBUG_FLAGS "-std=c++11 -g -O1 -L/usr/X11R6/lib -lm -lpthread -lX11 -pedantic")
	set(RELEASE_FLAGS "-std=c++11 -O3 -L/usr/X11R6/lib -lm -lpthread -lX11 -pedantic")
else()
	set(DEBUG_FLAGS "-std=c++11 -g -O1 -pedantic")
	set(RELEASE_FLAGS "-std=c++11 -O3 -pedantic")
endif()

set(CMAKE_CXX_FLAGS ${RELEASE_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})
set(CMAKE_CONFIGURATION_TYPES Debug Release)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()

add_subdirectory(test)