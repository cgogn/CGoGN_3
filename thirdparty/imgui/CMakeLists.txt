project(ImGUI LANGUAGES C CXX)

find_package(glfw3 REQUIRED)

set(HEADER_FILES
imconfig.h
imgui.h
imgui_impl_glfw.h
imgui_impl_opengl3.h
imgui_internal.h
imstb_rectpack.h
imstb_textedit.h
imstb_truetype.h
)

set(SOURCE_FILES
gl3w.c
imgui.cpp
imgui_draw.cpp
imgui_impl_glfw.cpp
imgui_impl_opengl3.cpp
imgui_widgets.cpp
)


add_library(${PROJECT_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cgogn/thirdparty>
)

target_link_libraries(${PROJECT_NAME} PUBLIC glfw)

set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")

install(FILES ${HEADER_FILES} DESTINATION "include/cgogn/thirdparty")
install(FILES ${HEADER_FILES} DESTINATION "include/cgogn/thirdparty")

#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/cgogn_rendering_export.h
#    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cgogn/rendering
#    COMPONENT cgogn_rendering_headers
#)



cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")
