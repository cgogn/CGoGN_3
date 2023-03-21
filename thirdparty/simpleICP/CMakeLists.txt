cmake_minimum_required(VERSION 3.15)
project(simpleICP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

set(HEADER_FILES
        simpleicp.h
        pointcloud.h
        corrpts.h
        nanoflann.hpp
)

set(SOURCE_FILES
        simpleicp.cpp
        pointcloud.cpp
        corrpts.cpp
)

add_library(${PROJECT_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cgogn/thirdparty>
)

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

if(MSVC)
        target_compile_options(${PROJECT_NAME} PRIVATE "/bigobj")
endif()

install(FILES ${HEADER_FILES} DESTINATION "include/cgogn/thirdparty")

cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER thirdparty)