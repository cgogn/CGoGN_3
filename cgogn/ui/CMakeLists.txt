project(cgogn_ui
	LANGUAGES CXX C
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_rendering REQUIRED)
find_package(ImGUI REQUIRED)
find_package(synapse REQUIRED)

add_library(${PROJECT_NAME} STATIC "")
add_library(cgogn::ui ALIAS ${PROJECT_NAME})

target_sources(${PROJECT_NAME}
	PRIVATE
		"${CMAKE_CURRENT_LIST_DIR}/app.h"
		"${CMAKE_CURRENT_LIST_DIR}/app.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/camera.h"
		"${CMAKE_CURRENT_LIST_DIR}/camera.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/gl_viewer.h"
		"${CMAKE_CURRENT_LIST_DIR}/gl_viewer.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/inputs.h"
		"${CMAKE_CURRENT_LIST_DIR}/module.h"
		"${CMAKE_CURRENT_LIST_DIR}/module.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/moving_frame.h"
		"${CMAKE_CURRENT_LIST_DIR}/view.h"
		"${CMAKE_CURRENT_LIST_DIR}/view.cpp"
)

set_target_properties(${PROJECT_NAME} PROPERTIES
	DEBUG_POSTFIX "_d"
	EXPORT_NAME ui
)

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}/thirdparty>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}/thirdparty/synapse/include>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cgogn/thirdparty>
)

target_link_libraries(${PROJECT_NAME} cgogn::core cgogn::rendering ImGUI synapse)

# Write out cgogn_ui_export.h to the current binary directory
generate_export_header(${PROJECT_NAME})

set(PKG_CONFIG_REQUIRES "cgogn_core cgogn_rendering ImGUI synapse")
configure_file(${PROJECT_SOURCE_DIR}/cgogn_ui.pc.in ${CMAKE_CURRENT_BINARY_DIR}/cgogn_ui.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_ui.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
	COMPONENT cgogn_ui_libraries
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/ui
	COMPONENT cgogn_ui_headers
	FILES_MATCHING PATTERN "*.h"
	REGEX "(examples|tests)" EXCLUDE
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_ui_export.h
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/ui
	COMPONENT cgogn_ui_headers
)

cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER cgogn)





add_subdirectory(modules/cmap_provider)
add_subdirectory(modules/surface_differential_properties)
add_subdirectory(modules/surface_filtering)
add_subdirectory(modules/surface_render)
