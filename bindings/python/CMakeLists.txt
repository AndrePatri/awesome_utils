cmake_minimum_required(VERSION 3.5)

set(PYBIND11_PYTHON_VERSION 3.8 CACHE STRING "")
find_package(pybind11 QUIET)

if(${pybind11_FOUND})
    
    message(STATUS "Python bindings will be compiled")
    include_directories(${CMAKE_CURRENT_SOURCE_DIR})
    pybind11_add_module(awesome_pyutils pyAwesomeUtils.cpp)
    set(PYBIND_LIB_LINKS awesome_utils)
    if(${WITH_XBOT2})
        list(APPEND PYBIND_LIB_LINKS
            xbot2::xbot2 xbot2::xbot2_ros_support)
    endif()

    target_link_libraries(awesome_pyutils PUBLIC ${PYBIND_LIB_LINKS})

    # get python install path
    execute_process(
        COMMAND python3 -c "if True:
            from distutils import sysconfig as sc
            print(sc.get_python_lib(prefix='', plat_specific=True))"
        OUTPUT_VARIABLE PYTHON_SITE
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Python install dir: ${PYTHON_SITE}")
    
    install(TARGETS awesome_pyutils
            DESTINATION ${PYTHON_SITE}/awesome_utils)
            
    file(WRITE ${CMAKE_BINARY_DIR}/__init__.py "")

    install(FILES  ${CMAKE_BINARY_DIR}/__init__.py DESTINATION ${PYTHON_SITE}/awesome_utils)
            
else()

    message(STATUS "Pybind11 not found, bindings cannot be compiled")

endif()
