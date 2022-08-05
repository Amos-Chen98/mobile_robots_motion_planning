cmake_minimum_required(VERSION 2.8.3)

project(lec5_hw)

set(CMAKE_CXX_FLAGS "-std=c++14")
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -fPIC")

find_package(Eigen3 REQUIRED)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  geometry_msgs
  visualization_msgs
)

include_directories(
    ${catkin_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    include
)

catkin_package()

add_executable(click_gen src/click_gen.cpp)

target_link_libraries(click_gen
  ${catkin_LIBRARIES}
)
