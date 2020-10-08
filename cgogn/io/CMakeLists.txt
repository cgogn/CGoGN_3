project(cgogn_io
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)

# Hide symbols by default
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)

find_package(Eigen3 3.2)

add_library(${PROJECT_NAME} STATIC "")
add_library(cgogn::io ALIAS ${PROJECT_NAME})

target_sources(${PROJECT_NAME}
	PRIVATE
		"${CMAKE_CURRENT_LIST_DIR}/graph/graph_import.h"
		"${CMAKE_CURRENT_LIST_DIR}/graph/graph_import.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/graph/cg.h"
		"${CMAKE_CURRENT_LIST_DIR}/graph/cgr.h"
		"${CMAKE_CURRENT_LIST_DIR}/graph/skel.h"

		"${CMAKE_CURRENT_LIST_DIR}/surface/surface_import.h"
		"${CMAKE_CURRENT_LIST_DIR}/surface/surface_import.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/surface/off.h"

		"${CMAKE_CURRENT_LIST_DIR}/volume/volume_import.h"
		"${CMAKE_CURRENT_LIST_DIR}/volume/volume_import.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/volume/mesh.h"
		"${CMAKE_CURRENT_LIST_DIR}/volume/meshb.h"
		"${CMAKE_CURRENT_LIST_DIR}/volume/tet.h"
		
		"${CMAKE_CURRENT_LIST_DIR}/utils.h"
)

set_target_properties(${PROJECT_NAME} PROPERTIES
	DEBUG_POSTFIX "_d"
	EXPORT_NAME io
	LINKER_LANGUAGE CXX
)

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cgogn/thirdparty>
)

# Write out cgogn_io_export.h to the current binary directory
generate_export_header(${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} cgogn::core Eigen3::Eigen Meshb)

set(PKG_CONFIG_REQUIRES "cgogn_core cgogn_geometry")
configure_file(${PROJECT_SOURCE_DIR}/cgogn_io.pc.in ${CMAKE_CURRENT_BINARY_DIR}/cgogn_io.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_io.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
	COMPONENT cgogn_io_libraries
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/io
	COMPONENT cgogn_io_headers
	FILES_MATCHING PATTERN "*.h"
	REGEX "(examples|tests)" EXCLUDE
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_io_export.h
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/io
	COMPONENT cgogn_io_headers
)

cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")

# VisualStudio nice source organisation
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER cgogn)
