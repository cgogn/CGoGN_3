cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(cgogn_rendering_examples
	LANGUAGES CXX
)

#find_package(cgogn_core REQUIRED)
#find_package(cgogn_ui REQUIRED)
#find_package(cgogn_io REQUIRED)
#find_package(cgogn_rendering REQUIRED)
#find_package(cgogn_simulation REQUIRED)

#set(CGOGN_TEST_PREFIX "test_")

#add_executable(shallow_water shallow_water.cpp)
#target_link_libraries(shallow_water
#	cgogn::core
#	cgogn::ui
#	cgogn::io
#	cgogn::rendering
#	cgogn::simulation
#	${CMAKE_DL_LIBS}
#)

#if(APPLE)
#	find_library(CORE_FOUNDATION CoreFoundation)
#	find_library(CARBON Carbon)
#	target_link_libraries(shallow_water
#		${CORE_FOUNDATION}
#		${CARBON}
#	)
#endif()
