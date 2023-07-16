cmake_minimum_required(VERSION 3.15)
project(Blossom5 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# find_package(Eigen3 3.3 REQUIRED NO_MODULE)

set(HEADER_FILES
        LCA.h
        PerfectMatching.h
        PMimplementation.h
        PQ.h
        timer.h
        MinCost/MinCost.h
        GEOM/GeomPerfectMatching.h
        GEOM/GPMkdtree.h
)

set(SOURCE_FILES
        misc.cpp
        PMduals.cpp
        PMexpand.cpp
        PMinit.cpp
        PMinterface.cpp
        PMmain.cpp
        PMrepair.cpp
        PMshrink.cpp
        MinCost/MinCost.cpp
        GEOM/GPMinit.cpp
        GEOM/GPMinterface.cpp
        GEOM/GPMkdtree.cpp
        GEOM/GPMmain.cpp
)

add_library(${PROJECT_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cgogn/thirdparty>
)

# target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

if(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE "/bigobj")
endif()

install(FILES ${HEADER_FILES} DESTINATION "include/cgogn/thirdparty")

cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER thirdparty)
