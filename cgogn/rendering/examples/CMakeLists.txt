cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(cgogn_rendering_examples
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)

find_package(cgogn_module_cmap_provider REQUIRED)
find_package(cgogn_module_surface_differential_properties REQUIRED)
find_package(cgogn_module_surface_render REQUIRED)

if(APPLE)
  find_library(CORE_FOUNDATION NAMES CoreFoundation)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -framework CoreFoundation")
endif()

set(CGOGN_TEST_PREFIX "test_")
set(CGOGN_TEST_MESHES_PATH "${CMAKE_SOURCE_DIR}/data/meshes/")
add_definitions("-DCGOGN_TEST_MESHES_PATH=${CGOGN_TEST_MESHES_PATH}")

add_executable(simple_viewer simple_viewer.cpp)
target_link_libraries(simple_viewer
	cgogn::core
	cgogn::ui
	cgogn::module::cmap_provider
	cgogn::module::surface_differential_properties
	cgogn::module::surface_render
)
