project(ConnectedComponnentITK)

cmake_minimum_required(VERSION 3.4)

find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

add_executable(ConnectedComponentITK ConnectedComponentITK.cxx)
add_executable(ShapeAttributes ShapeAttributes.cxx)

target_link_libraries(ConnectedComponentITK ${ITK_LIBRARIES})
target_link_libraries(ShapeAttributes ${ITK_LIBRARIES})
