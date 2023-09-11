project(cgogn_modeling
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_geometry REQUIRED)

# Hide symbols by default
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)

add_library(${PROJECT_NAME} STATIC "")
add_library(cgogn::modeling ALIAS ${PROJECT_NAME})

set(src_list
        "${CMAKE_CURRENT_LIST_DIR}/algos/decimation/decimation.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/decimation/cell_queue.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/decimation/edge_approximator.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/decimation/edge_queue_update.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/decimation/QEM_helper.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/convex_hull.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/graph_resampling.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/graph_resampling.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/algos/graph_to_hex.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/graph_to_hex.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/algos/incidenceGraph_to_hex.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/incidenceGraph_to_hex.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/algos/graph_utils.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/graph_utils.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/algos/mesh_repair.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/skeleton.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/subdivision.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/subdivision_utils.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/subdivision/basic.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/subdivision/surface_loop.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/subdivision/surface_catmull_clark.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/remeshing/blossom_quad_remeshing.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/remeshing/topstoc.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/remeshing/pliant_remeshing.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/volume_utils.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/volume_utils.cpp"

		"${CMAKE_CURRENT_LIST_DIR}/ui_modules/surface_modeling.h"
		"${CMAKE_CURRENT_LIST_DIR}/ui_modules/surface_deformation.h"
		"${CMAKE_CURRENT_LIST_DIR}/ui_modules/skeleton_extractor_module.h"
		"${CMAKE_CURRENT_LIST_DIR}/ui_modules/tubular_mesh_module.h"
		"${CMAKE_CURRENT_LIST_DIR}/ui_modules/volume_deformation.h"
		"${CMAKE_CURRENT_LIST_DIR}/ui_modules/volume_mr_modeling.h"
)
target_sources(${PROJECT_NAME} PRIVATE ${src_list})
source_group(TREE ${CMAKE_CURRENT_LIST_DIR} FILES ${src_list})

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

target_link_libraries(${PROJECT_NAME} cgogn::core cgogn::geometry)

# Write out cgogn_modeling_export.h to the current binary directory
generate_export_header(${PROJECT_NAME})

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

# VisualStudio nice source organisation
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER cgogn)
