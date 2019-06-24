project(cgogn_modeling
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)

# Hide symbols by default
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)

add_library(${PROJECT_NAME} SHARED "")
add_library(cgogn::modeling ALIAS ${PROJECT_NAME})

target_sources(${PROJECT_NAME}
	PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/algos/decimation/decimation.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/decimation/cell_queue.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/decimation/edge_approximator.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/decimation/edge_queue_edge_length.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/subdivision.h"
		"${CMAKE_CURRENT_LIST_DIR}/pipo.cpp"
)

if(${CGOGN_EXTERNAL_TEMPLATES})
	target_sources(${PROJECT_NAME}
		PRIVATE
			"${CMAKE_CURRENT_LIST_DIR}/external_templates.cpp"
	)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES
	DEBUG_POSTFIX "_d"
	EXPORT_NAME modeling
	LINKER_LANGUAGE CXX
)

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Write out cgogn_modeling_export.h to the current binary directory
generate_export_header(${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} cgogn::core cgogn::geometry)

set(PKG_CONFIG_REQUIRES "cgogn_core cgogn_geometry")
configure_file(${PROJECT_SOURCE_DIR}/cgogn_modeling.pc.in ${CMAKE_CURRENT_BINARY_DIR}/cgogn_modeling.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_modeling.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
	COMPONENT cgogn_modeling_libraries
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/modeling
	COMPONENT cgogn_modeling_headers
	FILES_MATCHING PATTERN "*.h"
	REGEX "(examples|tests)" EXCLUDE
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_modeling_export.h
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/modeling
	COMPONENT cgogn_modeling_headers
)


cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER cgogn)
