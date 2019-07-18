cmake_minimum_required(VERSION 3.7.2 FATAL_ERROR)

project(cgogn_geometry_examples
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)

find_package(cgogn_module_mesh_provider REQUIRED)
find_package(cgogn_module_surface_differential_properties REQUIRED)
find_package(cgogn_module_surface_filtering REQUIRED)
find_package(cgogn_module_surface_render REQUIRED)

set(CGOGN_TEST_PREFIX "test_")
set(CGOGN_TEST_MESHES_PATH "${CMAKE_SOURCE_DIR}/data/meshes/")
add_definitions("-DCGOGN_TEST_MESHES_PATH=${CGOGN_TEST_MESHES_PATH}")

add_executable(filtering filtering.cpp)
target_link_libraries(filtering
	cgogn::core
	cgogn::ui
	cgogn::module::mesh_provider
	cgogn::module::surface_differential_properties
	cgogn::module::surface_filtering
	cgogn::module::surface_render
)
