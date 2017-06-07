message(STATUS "Building tests")

set(target optics_tests)
message(STATUS "Building application ${target}")

set(sources
	test_main.cpp
)

add_executable(${target} ${sources})
add_test(all "${EXECUTABLE_OUTPUT_PATH}/${target}")