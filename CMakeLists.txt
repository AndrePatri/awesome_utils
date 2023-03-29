cmake_minimum_required(VERSION 3.5)

project(awesome_utils VERSION 1.0.0 LANGUAGES CXX)

option(WITH_XBOT2 "Compile and install xbot2-dependent components" FALSE)
if(${WITH_XBOT2})
    message(STATUS "XBot2-dependent components will be built")
    add_definitions(-DWITH_XBOT2) # can be used by preprocessor macros
endif()
option(WITH_MODEL_INTERFACE "Compile and install the model interface and all dependent components" FALSE)
if(${WITH_MODEL_INTERFACE})
    message(STATUS "The model interface and all dependent components will be built")
    add_definitions(-DWITH_MODEL_INTERFACE) # can be used by preprocessor macros
endif()
option(WITH_PYTHON "Compile and install python bindings" FALSE)
if(${WITH_PYTHON})
    message(STATUS "Python bindings will be built")
    add_definitions(-DWITH_PYTHON) # can be used by preprocessor macros
endif()
option(WITH_DOCS "Compile documentation" FALSE)
if(${WITH_DOCS})
    message(STATUS "Documentation will be built")
    add_definitions(-DWITH_DOCS) # can be used by preprocessor macros
endif()
option(BUILD_TESTS "Compile ${PROJECT_NAME} tests" FALSE)
if(${BUILD_TESTS})
    message(STATUS "Tests will be built")
    add_definitions(-DBUILD_TESTS) # can be used by preprocessor macros
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(Eigen3 REQUIRED)
find_package(matlogger2 REQUIRED)

find_package(pinocchio QUIET)

# base utilities
set(ORIENTATION_TRGT src/orientation_utils.cpp)
set(SIGN_PROC_UTILS_TRGT src/sign_proc_utils.cpp)
set(MODEL_INTERFACE_TRGT src/model_interface.cpp)

if(${WITH_XBOT2})

    set(XBOT2_UTILS_TRGT src/xbot2_utils.cpp)
    find_package(xbot_msgs REQUIRED)
    find_package(xbot2 REQUIRED)
    set(EC_CLIENT_SKIP_ECAT_MASTER TRUE) #hack to avoid find_package error on machines
    # where the ecat_master is not installed (most of them, since it's usually only on the
    # real robot)
    find_package(ec_xbot2_client QUIET)
    if(ec_xbot2_client_FOUND)
      add_definitions(-DEC_XBOT2_CLIENT_FOUND) # to be used by xbot2_utils executable
      # this option will expose methods which employ messages defined by xbot2_client
      # and read directly from XBot2 internal topics, for minimum latency.
      # To be improved: message definitions inside ec_xbot2_client should be put in a
      # separate package
    endif()

endif()

# more complex utilities (they employ one or more of the base utilities)
set(CALIB_UTILS_TRGT src/calib_utils.cpp)
set(TRAJ_UTILS_TRGT src/traj_utils.cpp)
set(CONTACT_EST_UTILS_TRGT src/contact_est_utils.cpp)

set(POWER_TRGT src/power_utils.cpp)
set(CARTESIAN_IMP_TRGT src/cartesian_imp_utils.cpp)


set(LIBRARY_TARGET_NAME awesome_utils) # set library name

set(${LIBRARY_TARGET_NAME}_SRC
    ${CALIB_UTILS_TRGT}
    ${SIGN_PROC_UTILS_TRGT}
    ${POWER_TRGT}
    ${TRAJ_UTILS_TRGT}
    ) # setting base utilities

set(${LIBRARY_TARGET_NAME}_INCLUDE_DIRECTORIES
    include/awesome_utils/
    include/xbot2_utils/
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${matlogger2_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}/)

set(${LIBRARY_TARGET_NAME}_INTERFACE_DIRECTORIES
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>"
    "${pinocchio_INCLUDE_DIRS}")

set(${LIBRARY_TARGET_NAME}_LIB_LINKS
    Eigen3::Eigen
    pinocchio::pinocchio
    matlogger2::matlogger2)

if(${WITH_XBOT2})# we add xbot2-dependent components

    list(APPEND ${LIBRARY_TARGET_NAME}_SRC
        ${XBOT2_UTILS_TRGT} )

    list(APPEND ${LIBRARY_TARGET_NAME}_INCLUDE_DIRECTORIES
        include/xbot2_utils/ ${xbot_msgs_INCLUDE_DIRS}/)

#    list(APPEND ${LIBRARY_TARGET_NAME}_INTERFACE_DIRECTORIES
#        include/xbot2_utils/ ${xbot_msgs_INCLUDE_DIRS}/)

    list(APPEND ${LIBRARY_TARGET_NAME}_LIB_LINKS
        xbot2::xbot2 xbot2::xbot2_ros_support)

endif()

if(${WITH_MODEL_INTERFACE}) # we add the model interface and all dependet-components

    list(APPEND ${LIBRARY_TARGET_NAME}_SRC
         ${MODEL_INTERFACE_TRGT} ${CONTACT_EST_UTILS_TRGT} ${CARTESIAN_IMP_TRGT} ${ORIENTATION_TRGT})

    list(APPEND ${LIBRARY_TARGET_NAME}_INCLUDE_DIRECTORIES
        ${pinocchio_INCLUDE_DIRS})

#    list(APPEND ${LIBRARY_TARGET_NAME}_INTERFACE_DIRECTORIES
#        include/xbot2_utils/ ${pinocchio_INCLUDE_DIRS}/)

    list(APPEND ${LIBRARY_TARGET_NAME}_LIB_LINKS
        pinocchio::pinocchio)

endif()

add_library(${LIBRARY_TARGET_NAME} SHARED ${${LIBRARY_TARGET_NAME}_SRC})

target_include_directories(${LIBRARY_TARGET_NAME} PRIVATE ${${LIBRARY_TARGET_NAME}_INCLUDE_DIRECTORIES})

target_link_libraries(${LIBRARY_TARGET_NAME} PRIVATE ${${LIBRARY_TARGET_NAME}_LIB_LINKS})
                                                     
target_compile_options(${LIBRARY_TARGET_NAME} PRIVATE -std=c++17)

set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES
        VERSION ${${PROJECT_NAME}_VERSION})

set(LIB_EXT ".so")
if(APPLE)
    set(LIB_EXT ".dylib")
endif()

target_compile_definitions(${LIBRARY_TARGET_NAME} PRIVATE -DAWESOMEUTILS_LIB_EXT="${LIB_EXT}")

target_include_directories(${LIBRARY_TARGET_NAME} INTERFACE ${${LIBRARY_TARGET_NAME}_INTERFACE_DIRECTORIES})

include_directories(${${LIBRARY_TARGET_NAME}_INCLUDE_DIRECTORIES}) # exposes include dirs also to downstream packages (i.e. py bindings and tests)

# Specify installation targets, typology and destination folders.
install(TARGETS  ${LIBRARY_TARGET_NAME}
    EXPORT   ${LIBRARY_TARGET_NAME}
    LIBRARY  DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
    ARCHIVE  DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
    RUNTIME  DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin)

# Mark cpp header files for installation
install(DIRECTORY include/${PROJECT_NAME}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    FILES_MATCHING PATTERN "*.h*"
    PATTERN ".svn" EXCLUDE)

if(${WITH_XBOT2})
    install(DIRECTORY include/"xbot2_utils"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
        FILES_MATCHING PATTERN "*.h*"
        PATTERN ".svn" EXCLUDE)
endif()

include(InstallBasicPackageFiles)
install_basic_package_files(${LIBRARY_TARGET_NAME}
                            VERSION ${${PROJECT_NAME}_VERSION}
                            COMPATIBILITY AnyNewerVersion
                            EXPORT ${LIBRARY_TARGET_NAME}
                            VARS_PREFIX ${LIBRARY_TARGET_NAME}
                            NO_CHECK_REQUIRED_COMPONENTS_MACRO)

# Compile python bindings
if(${WITH_PYTHON})
    add_subdirectory(bindings/python/)
endif()

if(${WITH_DOCS})
    add_subdirectory(doc/)
endif()

# Add the uninstall target
include(AddUninstallTarget)

# Add Testing target
if(BUILD_TESTS)
  enable_testing()
  add_custom_target(test_verbose ${CMAKE_CTEST_COMMAND} -V)
  add_subdirectory(tests)
endif()
