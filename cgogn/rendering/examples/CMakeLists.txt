cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(cgogn_rendering_examples
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering REQUIRED)

set(CGOGN_TEST_PREFIX "test_")
set(CGOGN_TEST_MESHES_PATH "${CMAKE_SOURCE_DIR}/data/meshes/")
add_definitions("-DCGOGN_TEST_MESHES_PATH=${CGOGN_TEST_MESHES_PATH}")

add_executable(simple_graph_viewer simple_graph_viewer.cpp)
target_link_libraries(simple_graph_viewer
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
)

add_executable(simple_surface_viewer simple_surface_viewer.cpp)
target_link_libraries(simple_surface_viewer
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
)

add_executable(simple_volume_viewer simple_volume_viewer.cpp)
target_link_libraries(simple_volume_viewer
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
)

if(APPLE)
	find_library(CORE_FOUNDATION CoreFoundation)
	find_library(CARBON Carbon)
	target_link_libraries(simple_graph_viewer 
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(simple_surface_viewer 
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(simple_volume_viewer 
		${CORE_FOUNDATION}
		${CARBON}
	)
endif()
