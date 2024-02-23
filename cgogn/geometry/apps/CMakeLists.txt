cmake_minimum_required(VERSION 3.7.2 FATAL_ERROR)

project(cgogn_geometry_examples
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering REQUIRED)

set(CGOGN_TEST_PREFIX "test_")


add_executable(filtering filtering.cpp)
target_link_libraries(filtering
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
)

add_executable(heat_method heat_method.cpp)
target_link_libraries(heat_method
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
)

add_executable(sphere_fitting sphere_fitting.cpp)
target_link_libraries(sphere_fitting
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
)

add_executable(registration registration.cpp)
target_link_libraries(registration
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	simpleICP
	${CMAKE_DL_LIBS}
)

add_executable(constriction constriction.cpp)
target_link_libraries(constriction
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	${CMAKE_DL_LIBS}
)

if(APPLE)
	find_library(CORE_FOUNDATION CoreFoundation)
	find_library(CARBON Carbon)
	target_link_libraries(filtering
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(heat_method
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(registration
		${CORE_FOUNDATION}
		${CARBON}
	)
endif()
