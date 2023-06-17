cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(cgogn_rendering_examples
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering REQUIRED)

set(CGOGN_TEST_PREFIX "test_")


if(APPLE)
	find_library(CORE_FOUNDATION CoreFoundation)
	find_library(CARBON Carbon)
	set (APPLE_LINK "${CORE_FOUNDATION} ${CARBON}")
endif()

add_executable(simple_gm2_viewer simple_gm2_viewer.cpp)
target_link_libraries(simple_gm2_viewer
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
	${APPLE_LINK}
)


add_executable(simple_gm3_viewer simple_gm3_viewer.cpp)
target_link_libraries(simple_gm3_viewer
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
	${APPLE_LINK}
)


add_executable(simple_graph_viewer simple_graph_viewer.cpp)
target_link_libraries(simple_graph_viewer
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
	${APPLE_LINK}
)

add_executable(simple_incidence_graph_viewer simple_incidence_graph_viewer.cpp)
target_link_libraries(simple_incidence_graph_viewer
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
	${APPLE_LINK}
)


add_executable(simple_surface_viewer simple_surface_viewer.cpp)
target_link_libraries(simple_surface_viewer
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
	${APPLE_LINK}
)

add_executable(simple_volume_viewer simple_volume_viewer.cpp)
target_link_libraries(simple_volume_viewer
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
	${APPLE_LINK}
)

#add_executable(surface_volume_viewer surface_volume_viewer.cpp)
#target_link_libraries(surface_volume_viewer
#	cgogn::core
#	cgogn::ui
#	cgogn::io
#	cgogn::rendering
#	${CMAKE_DL_LIBS}
#)

add_executable(test_shape test_shape.cpp)
target_link_libraries(test_shape
	cgogn::core
	cgogn::ui
	cgogn::rendering
	${CMAKE_DL_LIBS}
	${APPLE_LINK}
)

#if(APPLE)
#	find_library(CORE_FOUNDATION CoreFoundation)
#	find_library(CARBON Carbon)
#	target_link_libraries(simple_graph_viewer
#		${CORE_FOUNDATION}
#		${CARBON}
#	)
#	target_link_libraries(simple_incidence_graph_viewer
#		${CORE_FOUNDATION}
#		${CARBON}
#	)
#	target_link_libraries(simple_surface_viewer
#		${CORE_FOUNDATION}
#		${CARBON}
#	)
#	target_link_libraries(simple_volume_viewer
#		${CORE_FOUNDATION}
#		${CARBON}
#	)
#	target_link_libraries(surface_volume_viewer
#		${CORE_FOUNDATION}
#		${CARBON}
#	)
#endif()
