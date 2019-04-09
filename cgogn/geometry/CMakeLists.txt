project(cgogn_geometry
	LANGUAGES CXX
)

# Hide symbols by default
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)

find_package(Eigen3 3.2)

add_library(${PROJECT_NAME} SHARED "")
add_library(cgogn::geometry ALIAS cgogn_geometry)

target_sources(${PROJECT_NAME}
	PRIVATE
	    "${CMAKE_CURRENT_LIST_DIR}/types/vector_traits.h"

		"${CMAKE_CURRENT_LIST_DIR}/functions/normal.h"
		"${CMAKE_CURRENT_LIST_DIR}/functions/vector_ops.h"

		"${CMAKE_CURRENT_LIST_DIR}/algos/centroid.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/filtering.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/normal.h"
		"${CMAKE_CURRENT_LIST_DIR}/algos/subdivision.h"
)

if(${CGOGN_EXTERNAL_TEMPLATES})
    target_sources(${PROJECT_NAME}
        PRIVATE
            "${CMAKE_CURRENT_LIST_DIR}/external_templates.cpp"
    )
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES
	DEBUG_POSTFIX "_d"
	EXPORT_NAME geometry
	LINKER_LANGUAGE CXX
)

target_compile_definitions(${PROJECT_NAME} PRIVATE "-DCGOGN_GEOMETRY_DLL_EXPORT")
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC
	#$<BUILD_INTERFACE:${CGOGN_THIRDPARTY_EIGEN3_INCLUDE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cgogn/thirdparty>
)

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Write out cgogn_geometry_export.h to the current binary directory
include(GenerateExportHeader)
generate_export_header(cgogn_geometry)

target_link_libraries(${PROJECT_NAME} cgogn::core Eigen3::Eigen)

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
