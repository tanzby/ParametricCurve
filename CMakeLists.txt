cmake_minimum_required(VERSION 3.10)
project(ParametricCurve)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Eigen3 REQUIRED) 
add_library(curve
    src/cubic_spline.cpp
    src/cubic_spline.hpp
    src/bezier.hpp
    src/numerical.hpp
    src/curve.hpp
)
target_link_libraries(curve Eigen3::Eigen)
target_include_directories(curve PUBLIC src)

option(build_example "build example" ON)
message(STATUS "build example: " ${build_example})

if (build_example)
    message(STATUS "configuring example ...")
    add_subdirectory(example)
endif ()