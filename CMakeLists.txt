CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(p4)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/modules)
find_package(Eigen3 REQUIRED)

set(${PROJECT_NAME}_HDR_LIST graph_data_structures.hpp)
set(${PROJECT_NAME}_SRC_LIST main.cpp graph_data_structures.cpp)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_HDR_LIST} ${${PROJECT_NAME}_SRC_LIST})
