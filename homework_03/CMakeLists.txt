project(space_movements)

cmake_minimum_required(VERSION 3.20)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

SET(PROJECT_SOURCE space_movements)

add_executable(space_movements space_movements.cpp)
target_link_libraries(space_movements Eigen3::Eigen)
