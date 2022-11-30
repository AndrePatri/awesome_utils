cmake_minimum_required(VERSION 3.5)

project(awesome_utils VERSION 1.0.0 LANGUAGES CXX)

option(WITH_XBOT2 "Compile and install xbot2-dependent comonents" TRUE)
option(WITH_PYTHON "Compile and install python bindings" FALSE)
option(BUILD_TESTS "Compile ${PROJECT_NAME} tests" FALSE)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(catkin REQUIRED sensor_msgs geometry_msgs std_msgs message_generation)
find_package(Eigen3 REQUIRED)
find_package(pinocchio REQUIRED)
find_package(matlogger2 REQUIRED)

if(${WITH_XBOT2})
    message(STATUS "The xbot2-dependent components will be built")
    set(XBOT2_UTILS_TRGT src/xbot2_utils.cpp)
    find_package(xbot_msgs REQUIRED)
endif()

set(LIBRARY_TARGET_NAME awesome_utils) # set library name
set(CALIB_UTILS_TRGT src/calib_utils.cpp)
set(SIGN_PROC_UTILS_TRGT src/sign_proc_utils.cpp)
set(TRAJ_UTILS_TRGT src/traj_utils.cpp)
set(CONTACT_EST_UTILS_TRGT src/contact_est_utils.cpp)
set(MODEL_INTERFACE_TRGT src/model_interface.cpp)

set(${LIBRARY_TARGET_NAME}_SRC
        ${CALIB_UTILS_TRGT}
        ${SIGN_PROC_UTILS_TRGT}
        ${TRAJ_UTILS_TRGT}
        ${CONTACT_EST_UTILS_TRGT}
        ${MODEL_INTERFACE_TRGT}
        ${XBOT2_UTILS_TRGT}
    )

add_library(${LIBRARY_TARGET_NAME} SHARED ${${LIBRARY_TARGET_NAME}_SRC})

include_directories(
    include/
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${EIGEN3_INCLUDE_DIR}/
    ${catkin_INCLUDE_DIRS}/
    ${xbot_msgs_INCLUDE_DIRS}/
    ${CMAKE_BINARY_DIR}/generated/
    ${matlogger2_INCLUDE_DIRS}
)

target_link_libraries(${LIBRARY_TARGET_NAME} PRIVATE ${catkin_LIBRARIES}
                                                     Eigen3::Eigen
                                                     pinocchio::pinocchio
                                                     matlogger2::matlogger2)

target_compile_options(${LIBRARY_TARGET_NAME} PRIVATE -std=c++17)

set_target_properties(${LIBRARY_TARGET_NAME} PROPERTIES
        VERSION ${${PROJECT_NAME}_VERSION})

set(LIB_EXT ".so")
if(APPLE)
    set(LIB_EXT ".dylib")
endif()

target_compile_definitions(${LIBRARY_TARGET_NAME} PRIVATE -DAWESOMEUTILS_LIB_EXT="${LIB_EXT}")


target_include_directories(${LIBRARY_TARGET_NAME} INTERFACE "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
                                                            "$<INSTALL_INTERFACE:$<INSTALL_PREFIX>/${CMAKE_INSTALL_INCLUDEDIR}>"
                                                            ${EIGEN3_INCLUDE_DIR}
                                                            ${xbot_msgs_INCLUDE_DIRS}
                                                            ${pinocchio_INCLUDE_DIRS}
                                                            ${matlogger2_INCLUDE_DIRS}
)

# Specify installation targets, typology and destination folders.
install(TARGETS  ${LIBRARY_TARGET_NAME}
    EXPORT   ${LIBRARY_TARGET_NAME}
    LIBRARY  DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
    ARCHIVE  DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT lib
    RUNTIME  DESTINATION "${CMAKE_INSTALL_BINDIR}" COMPONENT bin)

# Mark cpp header files for installation
install(DIRECTORY include/${PROJECT_NAME}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    FILES_MATCHING PATTERN "*.h*"
    PATTERN ".svn" EXCLUDE
)

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

# Add the uninstall target
include(AddUninstallTarget)

# Add Testing target
if(BUILD_TESTS)
  enable_testing()
  add_custom_target(test_verbose ${CMAKE_CTEST_COMMAND} -V)
  add_subdirectory(tests)
endif()
