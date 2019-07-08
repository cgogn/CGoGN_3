cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(cgogn_rendering_examples
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_geometry REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering REQUIRED)
find_package(cgogn_modeling REQUIRED)

if(APPLE)
  find_library(CORE_FOUNDATION NAMES CoreFoundation)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -framework CoreFoundation")
endif()

set(CGOGN_TEST_PREFIX "test_")
set(CGOGN_TEST_MESHES_PATH "${CMAKE_SOURCE_DIR}/data/meshes/")
add_definitions("-DCGOGN_TEST_MESHES_PATH=${CGOGN_TEST_MESHES_PATH}")

add_executable(simple_viewer simple_viewer.cpp)
target_link_libraries(simple_viewer cgogn::core cgogn::io cgogn::rendering)

#add_executable(viewer_topo viewer_topo.cpp)
#target_link_libraries(viewer_topo cgogn::core cgogn::io cgogn::rendering cgogn::modeling)

#add_executable(drawing drawing.cpp)
#target_link_libraries(drawing cgogn::rendering)

#add_executable(text_drawing text_drawing.cpp)
#target_link_libraries(text_drawing cgogn::rendering)

#add_executable(viewer_topo3 viewer_topo3.cpp)
#target_link_libraries(viewer_topo3 cgogn::core cgogn::io cgogn::rendering)
