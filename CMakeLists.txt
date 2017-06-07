message(STATUS "CMake started")

cmake_minimum_required(VERSION 2.8)

project(OPTICS-Clustering)

if (CMAKE_COMPILER_IS_GNUCXX)
	set(DEBUG_FLAGS "-std=c++11 -g -O1 -lX11 -pthread -Wall -Wextra -Werror -pedantic")
	set(RELEASE_FLAGS "-std=c++11 -O3 -lX11 -pthread -Wall -Wextra -Werror -pedantic")
else()
	set(DEBUG_FLAGS "-std=c++11 -g -O1 -Wall -Wextra -Werror -pedantic")
	set(RELEASE_FLAGS "-std=c++11 -O3 -Wall -Wextra -Werror -pedantic")
endif ()


set(CMAKE_CXX_FLAGS ${RELEASE_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})
set(CMAKE_CONFIGURATION_TYPES Debug Release)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()

add_subdirectory(test)