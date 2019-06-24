cmake_minimum_required(VERSION 3.7.2 FATAL_ERROR)

project(cgogn_modeling_examples
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering REQUIRED)
find_package(cgogn_rendering_pureGL REQUIRED)

set(CGOGN_TEST_PREFIX "test_")
set(CGOGN_TEST_MESHES_PATH "${CMAKE_SOURCE_DIR}/data/meshes/")
add_definitions("-DCGOGN_TEST_MESHES_PATH=${CGOGN_TEST_MESHES_PATH}")

if (CGOGN_USE_QT)
add_executable(subdivision subdivision.cpp)
target_link_libraries(subdivision cgogn::core cgogn::io cgogn::rendering)
endif()

add_executable(subdivision_imgui subdivision_imgui.cpp)
target_link_libraries(subdivision_imgui cgogn::core cgogn::io cgogn::renderingPGL)

set_target_properties(subdivision PROPERTIES FOLDER examples/modeling)
