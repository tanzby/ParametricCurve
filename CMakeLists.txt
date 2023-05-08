cmake_minimum_required(VERSION 3.10)
project(ParametricCurve)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

include_directories(include)

option(build_example "build example" ON)
message(STATUS "build example: " ${build_example})

if (build_example)
    message(STATUS "configuring example ...")
    add_subdirectory(example)
endif ()