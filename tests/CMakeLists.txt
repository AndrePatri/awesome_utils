cmake_minimum_required(VERSION 3.5)
include(ExternalProject)

project(awesome_utils_tests LANGUAGES C CXX VERSION 1.0.0)

add_compile_options("-std=c++17")

# compite GTest
find_package(GTest)
find_package(Threads)

set(URDF_PATH ${CMAKE_CURRENT_SOURCE_DIR}/urdf/)

set(CMAKE_DEF_CONFIG_FNAME cmake_config) # set config file name
set(CMAKE_DEF_CONFIG_GENDIR ${CMAKE_CURRENT_BINARY_DIR}/generated) # set config file name
configure_file(cmake/${CMAKE_DEF_CONFIG_FNAME}.h.in ${CMAKE_DEF_CONFIG_GENDIR}/${CMAKE_DEF_CONFIG_FNAME}.h)

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

add_executable(TestModelInterface test_model_interface.cpp)
target_link_libraries(TestModelInterface ${TestLibs})
add_test(NAME TestModelInterface COMMAND TestModelInterface)

add_executable(TestContactEst test_contact_est.cpp)
target_link_libraries(TestContactEst ${TestLibs})
add_test(NAME TestContactEst COMMAND TestContactEst)

add_executable(TestCartImp test_cart_imp_utils.cpp)
target_link_libraries(TestCartImp ${TestLibs})
add_test(NAME TestCartImp COMMAND TestCartImp)

install(FILES ${CMAKE_DEF_CONFIG_GENDIR}/${CMAKE_DEF_CONFIG_FNAME}.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)

#*** python bindings tests ***#

if(${WITH_PYTHON})

    # retrieve python install path
    execute_process(
        COMMAND python -c "if True:
            from distutils import sysconfig as sc
            print(sc.get_python_lib(prefix='', plat_specific=True))"
        OUTPUT_VARIABLE PYTHON_SITE
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "******************** *Python install dir: ${PYTHON_SITE}")

    #install python bindings test files
    install(FILES pybind_tests/test_model_int_pybinds.py DESTINATION ${PYTHON_SITE}/awesome_utils)

endif()



