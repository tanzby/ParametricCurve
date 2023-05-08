find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

add_executable(mini_bezier_test mini_bezier_test.cpp)
target_link_libraries(mini_bezier_test curve)

add_executable(arc_length_test arc_length_test.cpp)
target_link_libraries(arc_length_test curve)

add_executable(fit_bezier_curve fit_bezier_curve.cpp)
target_link_libraries(fit_bezier_curve ${OpenCV_LIBS} curve)

add_executable(fit_cubic_spline fit_cubic_spline.cpp)
target_link_libraries(fit_cubic_spline ${OpenCV_LIBS} curve)

add_executable(arc_length_parameterized_bezier arc_length_parameterized_bezier.cpp)
target_link_libraries(arc_length_parameterized_bezier ${OpenCV_LIBS} curve)

add_executable(simple_cubic_interpolation simple_cubic_interpolation.cpp)
target_link_libraries(simple_cubic_interpolation ${OpenCV_LIBS} curve)
