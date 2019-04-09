
# the minimum version allowed for doxygen
set(DOXYGEN_MIN_VERSION 1.8.0)

# check doxygen
find_package(Doxygen ${DOXYGEN_MIN_VERSION} QUIET)
if(NOT DOXYGEN_FOUND)
	message(WARNING "Doxygen not found")
	unset(DOXYGEN_EXECUTABLE CACHE)
	return()
endif()

#target for generating a full documentation 

set(DOC_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR})

configure_file(full.dox full.dox @ONLY)

add_custom_target(
	doc-full
	COMMAND ${CMAKE_COMMAND} -E remove_directory ${DOC_OUTPUT_DIR}/full
	COMMAND ${DOXYGEN_EXECUTABLE} full.dox
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
	COMMENT "Generate full API documentation for ${PROJECT_NAME}"
	)

set_property(
	DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${DOC_OUTPUT_DIR}/full
	)

#install documentation
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/full/html DESTINATION doc/full COMPONENT doc-full OPTIONAL)
