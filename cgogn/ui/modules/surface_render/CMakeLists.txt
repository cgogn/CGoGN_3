project(cgogn_module_surface_render
	LANGUAGES CXX C
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)
find_package(cgogn_rendering REQUIRED)
find_package(cgogn_module_cmap_provider REQUIRED)
find_package(ImGUI REQUIRED)

add_library(${PROJECT_NAME} STATIC "")
add_library(cgogn::module::surface_render ALIAS ${PROJECT_NAME})

target_sources(${PROJECT_NAME}
	PRIVATE
		"${CMAKE_CURRENT_LIST_DIR}/surface_render.h"
		"${CMAKE_CURRENT_LIST_DIR}/surface_render.cpp"
)

set_target_properties(${PROJECT_NAME} PROPERTIES
	DEBUG_POSTFIX "_d"
	EXPORT_NAME module_surface_render
)

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}/thirdparty>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cgogn/thirdparty>
)

target_link_libraries(${PROJECT_NAME}
	cgogn::core
	cgogn::ui
	cgogn::rendering
	cgogn::module::cmap_provider
	ImGUI
)

# Write out cgogn_module_surface_render_export.h to the current binary directory
generate_export_header(${PROJECT_NAME})

set(PKG_CONFIG_REQUIRES "cgogn_core cgogn_rendering cgogn_module_cmap_provider ImGUI")
configure_file(${PROJECT_SOURCE_DIR}/cgogn_module_surface_render.pc.in ${CMAKE_CURRENT_BINARY_DIR}/cgogn_module_surface_render.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_module_surface_render.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
	COMPONENT cgogn_module_surface_render_libraries
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/module_surface_render
	COMPONENT cgogn_module_surface_render_headers
	FILES_MATCHING PATTERN "*.h"
	REGEX "(examples|tests)" EXCLUDE
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_module_surface_render_export.h
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/module_surface_render
	COMPONENT cgogn_ui_headers
)

cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER cgogn)
