project(cgogn_rendering_pureGL
	LANGUAGES CXX C
)

add_library(${PROJECT_NAME} SHARED "")

add_library(cgogn::renderingPGL ALIAS cgogn_rendering_pureGL)

target_sources(${PROJECT_NAME}
	PRIVATE
#	    "${CGOGN_THIRDPARTY_DIR}/imgui/gl3w.c"
        "${CMAKE_CURRENT_LIST_DIR}/drawer.h"
		"${CMAKE_CURRENT_LIST_DIR}/drawer.cpp"
#		"${CMAKE_CURRENT_LIST_DIR}/topo_drawer.h"
#		"${CMAKE_CURRENT_LIST_DIR}/topo_drawer.cpp"
#		"${CMAKE_CURRENT_LIST_DIR}/volume_drawer.h"
#		"${CMAKE_CURRENT_LIST_DIR}/volume_drawer.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/mesh_render.h"
		"${CMAKE_CURRENT_LIST_DIR}/mesh_render.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/wall_paper.h"
		"${CMAKE_CURRENT_LIST_DIR}/wall_paper.cpp"
#		"${CMAKE_CURRENT_LIST_DIR}/frame_manipulator.h"
#		"${CMAKE_CURRENT_LIST_DIR}/frame_manipulator.cpp"
#		"${CMAKE_CURRENT_LIST_DIR}/text_drawer.h"
#		"${CMAKE_CURRENT_LIST_DIR}/text_drawer.cpp"
#		"${CMAKE_CURRENT_LIST_DIR}/fonte.qrc"
        "${CMAKE_CURRENT_LIST_DIR}/types.h"
		"${CMAKE_CURRENT_LIST_DIR}/gl_image.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/vbo.h"
		"${CMAKE_CURRENT_LIST_DIR}/ebo.h"
		"${CMAKE_CURRENT_LIST_DIR}/vao.h"
		"${CMAKE_CURRENT_LIST_DIR}/texture.h"
		"${CMAKE_CURRENT_LIST_DIR}/fbo.h"
		"${CMAKE_CURRENT_LIST_DIR}/fbo.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/vbo_update.h"
		"${CMAKE_CURRENT_LIST_DIR}/mframe.h"
		"${CMAKE_CURRENT_LIST_DIR}/camera.h"
		"${CMAKE_CURRENT_LIST_DIR}/camera.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/pure_gl_viewer.h"
		"${CMAKE_CURRENT_LIST_DIR}/pure_gl_viewer.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/imgui_viewer.h"
		"${CMAKE_CURRENT_LIST_DIR}/imgui_viewer.cpp"

		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_program.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_program.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_simple_color.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_simple_color.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_flat.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_flat.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_flat_cpv.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_flat_cpv.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_color_per_vertex.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_color_per_vertex.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_scalar_per_vertex.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_scalar_per_vertex.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_vector_per_vertex.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_vector_per_vertex.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_phong.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_phong.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_phong_color.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_phong_color.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_bold_line.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_bold_line.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_bold_line_color.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_bold_line_color.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_round_point.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_round_point.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_round_point_color.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_round_point_color.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_point_sprite.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_point_sprite.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_explode_volumes.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_explode_volumes.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_explode_volumes_color.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_explode_volumes_color.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_explode_volumes_line.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_explode_volumes_line.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_texture.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_texture.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_text.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_text.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/shaders/shader_fullscreen_texture.h"
        "${CMAKE_CURRENT_LIST_DIR}/shaders/shader_fullscreen_texture.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_pick.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_pick.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_frame2d.h"
		"${CMAKE_CURRENT_LIST_DIR}/shaders/shader_frame2d.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/shaders/shader_histo.h"
        "${CMAKE_CURRENT_LIST_DIR}/shaders/shader_histo.cpp"
)

#[[
if (Qt5Gui_VERSION VERSION_GREATER 5.5.99)
	target_sources(${PROJECT_NAME}
		PRIVATE
		"${CMAKE_CURRENT_LIST_DIR}/transparency_shaders/shader_transparent_quad.h"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_shaders/shader_transparent_quad.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_shaders/shader_transparent_flat.h"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_shaders/shader_transparent_flat.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_shaders/shader_transparent_phong.h"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_shaders/shader_transparent_phong.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_shaders/shader_transparent_volumes.h"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_shaders/shader_transparent_volumes.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_shaders/shader_copy_depth.h"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_shaders/shader_copy_depth.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_drawer.h"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_drawer.cpp"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_volume_drawer.h"
		"${CMAKE_CURRENT_LIST_DIR}/transparency_volume_drawer.cpp"
	)
else()
	message(WARNING "Qt version < 5.6 : not compiling transparency rendering (need MRT FBO)")
endif()
]]
set_target_properties(${PROJECT_NAME} PROPERTIES
										DEBUG_POSTFIX "_d"
										EXPORT_NAME rendering)

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}/thirdparty>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}/thirdparty/imgui>
	$<BUILD_INTERFACE:${CGOGN_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cgogn/thirdparty>
)

#target_link_libraries(${PROJECT_NAME} cgogn::core cgogn::geometry ImGUI)
target_link_libraries(${PROJECT_NAME} cgogn::core ImGUI)


generate_export_header(cgogn_rendering_pureGL)

set(PKG_CONFIG_REQUIRES "cgogn_core cgogn_geometry")
configure_file(${PROJECT_SOURCE_DIR}/cgogn_rendering_pureGL.pc.in ${CMAKE_CURRENT_BINARY_DIR}/cgogn_rendering_pureGL.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_rendering_pureGL.pc
	DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
	COMPONENT ccgogn_rendering_pureGL_libraries
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/rendering_pureGL
	COMPONENT cgogn_rendering_pureGL_headers
	FILES_MATCHING PATTERN "*.h"
	REGEX "(examples|tests)" EXCLUDE
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_rendering_pureGL_export.h
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/rendering_pureGL
	COMPONENT cgogn_rendering_pureGL_headers
)

cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")



set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER cgogn)
