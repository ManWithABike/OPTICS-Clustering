
message(STATUS "Building Unit Tests ${UNITTEST}")

find_package(GTest REQUIRED)
find_package(GMock REQUIRED)

macro(_add_test _NAME)
    add_executable(${_NAME} ${_NAME}.cpp)
    add_test(NAME ${_NAME} COMMAND ${_NAME})

    target_include_directories(${_NAME} SYSTEM PUBLIC
                                ${GTest_INCLUDE_DIR}
                                ${GMock_INCLUDE_DIR}
                                )
    target_link_libraries(${_NAME} ${GTest_LIBRARIES}
                                    ${GMock_LIBRARIES}
                                    ${CMAKE_THREAD_LIBS_INIT}
                                    )
endmacro()

_add_test(reachability_test)


add_custom_target(unittest reachability_test
                        COMMENT "Running unittests\n\n"
                        VERBATIM
                        )

