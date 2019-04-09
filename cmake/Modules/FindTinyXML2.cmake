# FindTinyXML2.cmake
#
# Finds the TinyXML2 library
#
# This will define the following variables
#
#    	TinyXML2_FOUND - system has TinyXML2
#    	TinyXML2_INCLUDE_DIRS - the TinyXML2 include directory
#	 	TinyXML2_LIBRARIES - The libraries needed to use TinyXML2
#		TinyXML2_DEFINITIONS - Compiler switches required for using TinyXML2
#
# and the following imported target
#
#     TinyXML2::TinyXML2
#
# Author: Lionel Untereiner - lionel.untereiner@geosiris.com

find_package(PkgConfig)
pkg_check_modules(PC_TinyXML2 QUIET TinyXML2)
set(TinyXML2_DEFINITIONS ${PC_TinyXML2_CFLAGS_OTHER})

find_path(TinyXML2_INCLUDE_DIR
    NAMES tinyxml2.h
    HINTS ${PC_TinyXML2_INCLUDEDIR} ${PC_TinyXML2_INCLUDE_DIRS}
	DOC "The TinyXML2 include directory"
)

if(NOT TinyXML2_LIBRARY)
	find_library(TinyXML2_LIBRARY_RELEASE
		NAMES tinyxml2
		HINTS ${PC_TinyXML2_LIBDIR} ${PC_TinyXML2_LIBRARY_DIRS}
		DOC "The release TinyXML2 library"
	)

	find_library(TinyXML2_LIBRARY_DEBUG
		NAMES tinyxml2d
		HINTS ${PC_TinyXML2_LIBDIR} ${PC_TinyXML2_LIBRARY_DIRS}
		DOC "The debug TinyXML2 library"
	)
endif()

include(SelectLibraryConfigurations)
select_library_configurations(TinyXML2)

set(TinyXML2_VERSION ${PC_TinyXML2_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TinyXML2
    REQUIRED_VARS TinyXML2_INCLUDE_DIR TinyXML2_LIBRARY
    VERSION_VAR TinyXML2_VERSION
)

if(TinyXML2_FOUND)
    set(TinyXML2_INCLUDE_DIRS ${TinyXML2_INCLUDE_DIR})
	set(TinyXML2_LIBRARIES ${TinyXML2_LIBRARY})

	if(NOT TARGET TinyXML2::TinyXML2)
	    add_library(TinyXML2::TinyXML2 UNKNOWN IMPORTED)
	    set_target_properties(TinyXML2::TinyXML2 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${TinyXML2_INCLUDE_DIR}")

	    if(TinyXML2_LIBRARY_RELEASE)
        	set_property(TARGET TinyXML2::TinyXML2 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
        	set_target_properties(TinyXML2::TinyXML2 PROPERTIES IMPORTED_LOCATION_RELEASE "${TinyXML2_LIBRARY_RELEASE}")
		endif()

		if(TinyXML2_LIBRARY_DEBUG)
        	set_property(TARGET TinyXML2::TinyXML2 APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
        	set_target_properties(TinyXML2::TinyXML2 PROPERTIES IMPORTED_LOCATION_DEBUG "${TinyXML2_LIBRARY_DEBUG}")
		endif()

		if(NOT TinyXML2_LIBRARY_RELEASE AND NOT TinyXML2_LIBRARY_DEBUG)
        	set_property(TARGET TinyXML2::TinyXML2 APPEND PROPERTY IMPORTED_LOCATION "${TinyXML2_LIBRARY}")
		endif()
	endif()
endif()

mark_as_advanced(TinyXML2_INCLUDE_DIR TinyXML2_LIBRARY)
