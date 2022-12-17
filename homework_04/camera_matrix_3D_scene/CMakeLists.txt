project(camera)

cmake_minimum_required(VERSION 3.20)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(OpenCV REQUIRED)

SET(PROJECT_SOURCE main.cpp camera.cpp)
SET(PROJECT_HEADER camera.h)

add_executable(camera main.cpp camera.cpp)
target_link_libraries(camera Eigen3::Eigen ${OpenCV_LIBS})