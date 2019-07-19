
include(CMakeParseArguments)

#!
# @brief automatic deduction of CMAKE_BUILD_TYPE from CMAKE_CURRENT_BINARY_DIR
# @details
# if CMAKE_CURRENT_BINARY_DIR end with Debug/debug/DEBUG, set Debug mode else Release mode.
#
function(deduce_build_type)
set(BIN_DIR "")
string(TOLOWER ${CMAKE_CURRENT_BINARY_DIR} BIN_DIR)
if (${BIN_DIR} MATCHES "(.*)debug")
	set(CMAKE_BUILD_TYPE "Debug" PARENT_SCOPE)
else()
	if (${BIN_DIR} MATCHES "(.*)release")
		set(CMAKE_BUILD_TYPE "Release" PARENT_SCOPE)
	endif()
endif()
endfunction(deduce_build_type)

##############################################################################################
#                                 cgogn_create_package macro                                 #
# This macro is a helper to create package configuration and version files. These files are  #
# needed when using the find_package command.                                                #
# This macro generate 2 versions of each file : one for the build tree and another for the   #
# install tree.                                                                              #
# Build tree:                                                                                #
# 1.<build-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>Targets.cmake             #
# 2.<build-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>Config.cmake              #
# 3.<build-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>ConfigVersion.cmake       #
#                                                                                            #
# Install tree:                                                                              #
# 1.<install-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>Targets.cmake           #
# 2.<install-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>Config.cmake            #
# 3.<install-dir>/lib/cmake/<cmake-project-name>/<cmake-project-name>ConfigVersion.cmake     #
#                                                                                            #
# Usage example : find_package(cgogn_core); find_package(cgogn_io)                           #
# Note: template config files are located in  the cmake/ConfigFiles directory.               #
# By convention they have to define the following two variables:                             #
# cmake/<cmake-project-name>_LIBRARIES                                                       #
# cmake/<cmake-project-name>_INCLUDE_DIRS                                                    #
##############################################################################################

macro(cgogn_create_package package_root_dir)

######## 1. Build tree

export(TARGETS ${PROJECT_NAME}
	NAMESPACE cgogn::
	FILE "${CMAKE_BINARY_DIR}/lib/cmake/${PROJECT_NAME}/${PROJECT_NAME}Targets.cmake"
)

export(PACKAGE ${PROJECT_NAME})

configure_package_config_file(
	"${package_root_dir}/${PROJECT_NAME}Config.cmake.in"
	"${CMAKE_BINARY_DIR}/lib/cmake/${PROJECT_NAME}/${PROJECT_NAME}Config.cmake"
	INSTALL_DESTINATION "${CMAKE_BINARY_DIR}/lib/cmake/${PROJECT_NAME}"
	NO_SET_AND_CHECK_MACRO
	NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

write_basic_package_version_file(
	"${CMAKE_BINARY_DIR}/lib/cmake/${PROJECT_NAME}/${PROJECT_NAME}ConfigVersion.cmake"
	VERSION ${CGOGN_VERSION_MAJOR}.${CGOGN_VERSION_MINOR}.${CGOGN_VERSION_PATCH}
	COMPATIBILITY ExactVersion
)

######## 2. Install tree

install(TARGETS ${PROJECT_NAME}
	EXPORT ${PROJECT_NAME}Targets
	RUNTIME
		DESTINATION ${CMAKE_INSTALL_BINDIR}
		COMPONENT ${PROJECT_NAME}_applications
	LIBRARY
		DESTINATION ${CMAKE_INSTALL_LIBDIR}
		COMPONENT ${PROJECT_NAME}_libraries
	ARCHIVE
		DESTINATION ${CMAKE_INSTALL_LIBDIR}
		COMPONENT ${PROJECT_NAME}_libraries
)

install(EXPORT ${PROJECT_NAME}Targets
		NAMESPACE cgogn::
		DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
		COMPONENT ${PROJECT_NAME}_libraries
)

configure_package_config_file(
	"${package_root_dir}/${PROJECT_NAME}Config.cmake.in"
	"${CMAKE_BINARY_DIR}/cmake/${PROJECT_NAME}/${PROJECT_NAME}Config.cmake"
	INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
	NO_SET_AND_CHECK_MACRO
	NO_CHECK_REQUIRED_COMPONENTS_MACRO
)

install(FILES "${CMAKE_BINARY_DIR}/cmake/${PROJECT_NAME}/${PROJECT_NAME}Config.cmake"
	DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
	COMPONENT ${PROJECT_NAME}_libraries
)

write_basic_package_version_file(
	"${CMAKE_BINARY_DIR}/cmake/${PROJECT_NAME}/${PROJECT_NAME}ConfigVersion.cmake"
	VERSION ${CGOGN_VERSION_MAJOR}.${CGOGN_VERSION_MINOR}.${CGOGN_VERSION_PATCH}
	COMPATIBILITY ExactVersion
)

install(FILES "${CMAKE_BINARY_DIR}/cmake/${PROJECT_NAME}/${PROJECT_NAME}ConfigVersion.cmake"
	DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
	COMPONENT ${PROJECT_NAME}_libraries
)

endmacro()


##############################################################################################

function(cgogn_list_subdirectory result current_directory)
	file(GLOB children RELATIVE ${current_directory} ${current_directory}/*)
	set(dirlist "")
	foreach(child ${children})
		if(IS_DIRECTORY ${current_directory}/${child})
			list(APPEND dirlist ${child})
		endif()
	endforeach()
	set(${result} ${dirlist} PARENT_SCOPE)
endfunction()
