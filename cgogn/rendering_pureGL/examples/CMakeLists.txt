cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(cgogn_rendering_puregl_examples
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_geometry REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering_pureGL REQUIRED)
find_package(cgogn_modeling REQUIRED)

if(MACOSX)
  find_library(CORE_FOUNDATION NAMES CoreFoundation)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -framework CoreFoundation")
endif()

set(CGOGN_TEST_PREFIX "test_")
set(CGOGN_TEST_MESHES_PATH "${CMAKE_SOURCE_DIR}/data/meshes/")
add_definitions("-DCGOGN_TEST_MESHES_PATH=${CGOGN_TEST_MESHES_PATH}")

add_executable(simple_viewerPGL simple_viewer.cpp)
target_link_libraries(simple_viewerPGL cgogn::core cgogn::io cgogn::renderingPGL)


#add_executable(viewer_topoPGL viewer_topo.cpp)
#target_link_libraries(viewer_topoPGL cgogn::core cgogn::io cgogn::renderingPGL cgogn::modeling)

#add_executable(drawingPGL drawing.cpp)
#target_link_libraries(drawingPGL cgogn::renderingPGL)

#add_executable(text_drawingPGL text_drawing.cpp)
#target_link_libraries(text_drawingPGL cgogn::renderingPGL)

#add_executable(viewer_topo3PGL viewer_topo3.cpp)
#target_link_libraries(viewer_topo3PGL cgogn::core cgogn::io cgogn::renderingPGL)

#[[add_executable(viewer_per_face viewer_per_face.cpp)
target_link_libraries(viewer_per_face cgogn::core cgogn::io cgogn::rendering)


add_executable(picking_viewer picking_viewer.cpp)
target_link_libraries(picking_viewer cgogn::core cgogn::io cgogn::rendering)

add_executable(mt_viewer mt_viewer.cpp)
target_link_libraries(mt_viewer cgogn::core cgogn::io cgogn::rendering)

add_executable(point_set_drawing point_set_drawing.cpp)
target_link_libraries(point_set_drawing cgogn::rendering cgogn::io QOGLViewer)

add_executable(polyline_drawing polyline_drawing.cpp)
target_link_libraries(polyline_drawing cgogn::rendering cgogn::io QOGLViewer)


set_target_properties(simple_viewer viewer_per_face viewer_topo drawing picking_viewer viewer_topo3 mt_viewer text_drawing point_set_drawing
PROPERTIES FOLDER examples/rendering)


if (Qt5Gui_VERSION VERSION_GREATER 5.5.99)

	add_executable(transparency_viewer transparency_viewer.cpp)
	target_link_libraries(transparency_viewer cgogn::core cgogn::io cgogn::rendering)

	add_executable(transparencies_viewer transparencies_viewer.cpp)
	target_link_libraries(transparencies_viewer cgogn::core cgogn::io cgogn::rendering)

	add_executable(transparency_volume_viewer transparency_volume_viewer.cpp)
	target_link_libraries(transparency_volume_viewer cgogn::core cgogn::io cgogn::rendering)

	set_target_properties(transparency_viewer transparency_volume_viewer transparencies_viewer 
	PROPERTIES FOLDER examples/rendering)

endif()
]]
