cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(cgogn_rendering_examples
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering REQUIRED)
find_package(cgogn_modeling REQUIRED)

set(CGOGN_TEST_PREFIX "test_")

add_executable(modeling modeling.cpp)
target_link_libraries(modeling
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
)

add_executable(tubular_mesh tubular_mesh.cpp)
target_link_libraries(tubular_mesh
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
)

add_executable(deformation deformation.cpp)
target_link_libraries(deformation
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
)

add_executable(mr_volume mr_volume.cpp)
target_link_libraries(mr_volume
	cgogn::core
	cgogn::ui
	cgogn::io
	cgogn::rendering
	cgogn::modeling
	${CMAKE_DL_LIBS}
)

if(APPLE)
	find_library(CORE_FOUNDATION CoreFoundation)
	find_library(CARBON Carbon)
	target_link_libraries(modeling
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(tubular_mesh
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(deformation
		${CORE_FOUNDATION}
		${CARBON}
	)
	target_link_libraries(mr_volume
		${CORE_FOUNDATION}
		${CARBON}
	)
endif()
