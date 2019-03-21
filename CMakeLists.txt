cmake_minimum_required(VERSION 3.5.0)
project(Curve)

set(CMAKE_CXX_STANDARD 11)

find_package(PCL REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(example example.cpp)
target_link_libraries(example ${OpenCV_LIBS})