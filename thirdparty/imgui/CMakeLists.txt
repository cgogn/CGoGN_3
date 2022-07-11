project(ImGUI LANGUAGES C CXX)

find_package(glfw3 REQUIRED)

set(HEADER_FILES
	GL/gl3w.h
	GL/glcorearb.h
	imconfig.h
	imgui.h
	imgui_impl_glfw.h
	imgui_impl_opengl3.h
	imgui_impl_opengl3_loader.h
	imgui_internal.h
	imstb_rectpack.h
	imstb_textedit.h
	imstb_truetype.h
	imgui-knobs.h
)

set(SOURCE_FILES
	gl3w.c
	imgui.cpp
	imgui_demo.cpp
	imgui_draw.cpp
	imgui_tables.cpp
	imgui_impl_glfw.cpp
	imgui_impl_opengl3.cpp
	imgui_widgets.cpp
	imgui-knobs.cpp
)

add_library(${PROJECT_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cgogn/thirdparty>
)

target_compile_definitions(${PROJECT_NAME} PUBLIC IMGUI_IMPL_OPENGL_LOADER_GL3W NOMINMAX)
if (NOT MSVC)
    target_link_libraries(${PROJECT_NAME} PUBLIC glfw dl)
else()
    target_link_libraries(${PROJECT_NAME} PUBLIC glfw)
endif()

install(FILES ${HEADER_FILES} DESTINATION "include/cgogn/thirdparty")

cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")
