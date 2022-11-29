cmake_minimum_required(VERSION 3.5)

set(PYBIND11_PYTHON_VERSION 3.8 CACHE STRING "")
find_package(pybind11 QUIET)

add_compile_options(-std=c++14)

if(${pybind11_FOUND})
    
    message(STATUS "Will compile python bindings")
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    pybind11_add_module(awesomeutils pyAwesomeUtils.cpp)
    target_link_libraries(awesomeutils PUBLIC awesome_utils)
    
    # get python install path 
    execute_process(
        COMMAND python -c "if True:
            from distutils import sysconfig as sc
            print(sc.get_python_lib(prefix='', plat_specific=True))"
        OUTPUT_VARIABLE PYTHON_SITE
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Python install dir: ${PYTHON_SITE}")
    
    install(TARGETS awesomeutils
            DESTINATION ${PYTHON_SITE}/awesome_utils)
            
    file(WRITE ${CMAKE_BINARY_DIR}/__init__.py "")
    install(FILES  ${CMAKE_BINARY_DIR}/__init__.py DESTINATION ${PYTHON_SITE}/awesome_utils)

            
else()

    message(STATUS "Pybind11 not found, bindings won't be available")
    
endif()
