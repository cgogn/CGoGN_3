project(cgogn_geometry
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)

# Hide symbols by default
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)

find_package(Eigen3 3.2)

add_library(${PROJECT_NAME} STATIC "")
add_library(cgogn::geometry ALIAS ${PROJECT_NAME})

target_sources(${PROJECT_NAME}
	PRIVATE
	    "${CMAKE_CURRENT_LIST_DIR}/types/vector_traits.h"

		"${CMAKE_CURRENT_LIST_DIR}/functions/angle.h"
		"${CMAKE_CURRENT_LIST_DIR}/functions/area.h"
		"${CMAKE_CURRENT_LIST_DIR}/functions/inclusion.h"
		"${CMAKE_CURRENT_LIST_DIR}/functions/normal.h"
                "${CMAKE_CURRENT_LIST_DIR}/functions/distance.h"

        "${CMAKE_CURRENT_LIST_DIR}/algos/angle.h"
        "${CMAKE_CURRENT_LIST_DIR}/algos/area.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/centroid.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/curvature.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/filtering.h"
        "${CMAKE_CURRENT_LIST_DIR}/algos/length.h"
        "${CMAKE_CURRENT_LIST_DIR}/algos/normal.h"
        "${CMAKE_CURRENT_LIST_DIR}/algos/selection.h"
        "${CMAKE_CURRENT_LIST_DIR}/algos/ear_triangulation.h"

		"${CMAKE_CURRENT_LIST_DIR}/algos/pipo.cpp"
)

set_target_properties(${PROJECT_NAME} PROPERTIES
	DEBUG_POSTFIX "_d"
	EXPORT_NAME geometry
	LINKER_LANGUAGE CXX
)

target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cgogn/thirdparty>
)

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(${PROJECT_NAME} cgogn::core Eigen3::Eigen)

# Write out cgogn_geometry_export.h to the current binary directory
generate_export_header(${PROJECT_NAME})

set(PKG_CONFIG_REQUIRES "cgogn_core")
configure_file(${PROJECT_SOURCE_DIR}/cgogn_geometry.pc.in ${CMAKE_CURRENT_BINARY_DIR}/cgogn_geometry.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_geometry.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
	COMPONENT cgogn_geometry_libraries
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/geometry
	COMPONENT cgogn_geometry_headers
	FILES_MATCHING PATTERN "*.h"
	REGEX "(examples|tests)" EXCLUDE
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_geometry_export.h
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/geometry
	COMPONENT cgogn_geometry_headers
)

cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")

# VisualStudio nice source organisation
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER cgogn)
