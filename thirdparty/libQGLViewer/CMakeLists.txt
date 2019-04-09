set(CGOGN_THIRDPARTY_QOGLVIEWER_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "QOGLViewer include directory")

project(QOGLViewer)
find_package(OpenGL REQUIRED)
find_package(Qt5 5.4.0 COMPONENTS Widgets REQUIRED)

set(CMAKE_AUTOMOC ON)

# in order to be able to include generated files (ui and moc)
#include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(HEADER_FILES
	QOGLViewer/qoglviewer.h
	QOGLViewer/camera.h
	QOGLViewer/manipulatedFrame.h
	QOGLViewer/manipulatedCameraFrame.h
	QOGLViewer/frame.h
	QOGLViewer/keyFrameInterpolator.h
	QOGLViewer/config.h
	QOGLViewer/constraint.h
	QOGLViewer/mouseGrabber.h
	QOGLViewer/quaternion.h
	QOGLViewer/vec.h
)

set(SOURCE_FILES
	QOGLViewer/qoglviewer.cpp
	QOGLViewer/camera.cpp
	QOGLViewer/manipulatedFrame.cpp
	QOGLViewer/manipulatedCameraFrame.cpp
	QOGLViewer/frame.cpp
	QOGLViewer/saveSnapshot.cpp
	QOGLViewer/constraint.cpp
	QOGLViewer/keyFrameInterpolator.cpp
	QOGLViewer/mouseGrabber.cpp
	QOGLViewer/quaternion.cpp
	QOGLViewer/vec.cpp
)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC ${OPENGL_LIBRARY} Qt5::Widgets)

if(WIN32)
	add_definitions(-DCREATE_QGLVIEWER_DLL -DNOMINMAX /W3)
else()
	if(NOT MSVC)
		if (CGOGN_CPP_STD STREQUAL "11")
			target_compile_options(${PROJECT_NAME} PUBLIC "-std=c++11")
		elseif(CGOGN_CPP_STD STREQUAL "14")
			target_compile_options(${PROJECT_NAME} PUBLIC "-std=c++14")
		else(CGOGN_CPP_STD STREQUAL "11")
			target_compile_options(${PROJECT_NAME} PUBLIC "-std=c++17")
		endif(CGOGN_CPP_STD STREQUAL "11")
	endif(NOT MSVC)
endif()

# for glu warning on mac
target_compile_options(${PROJECT_NAME} PUBLIC
		$<$<CXX_COMPILER_ID:Clang>:-Wno-deprecated-declarations>)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/QOGLViewer"
	DESTINATION "include/cgogn/thirdparty"
	FILES_MATCHING PATTERN "*.h")

cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")
