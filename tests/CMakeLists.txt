cmake_minimum_required(VERSION 3.5)
include(ExternalProject)

project(awesome_utils_tests LANGUAGES C CXX VERSION 1.0.0)

add_compile_options("-std=c++20")

# compite GTest
find_package(GTest)
find_package(Threads)
if(NOT GTEST_FOUND)
    set(gTestSource "/usr/src/gtest")
    ExternalProject_Add(GTest-ext SOURCE_DIR ${gTestSource}
                                  PREFIX "${CMAKE_CURRENT_BINARY_DIR}/external"
                                  INSTALL_COMMAND "")
    set(GTEST_LIB_DIRS "${CMAKE_CURRENT_BINARY_DIR}/external/src/GTest-ext-build/")
    set(GTEST_INCLUDE_DIRS ${gTestSource})
    set(GTEST_BOTH_LIBRARIES gtest gtest_main)
endif()


set(TestLibs awesome_utils ${GTEST_BOTH_LIBRARIES} -lm -lpthread)
include_directories(${GTEST_INCLUDE_DIRS})
link_directories(${GTEST_LIB_DIRS})

add_executable(SmallTest small_test.cpp)
target_link_libraries(SmallTest ${TestLibs}))
add_dependencies(TestApi GTest-ext)
add_test(NAME SmallTest COMMAND SmallTest)

