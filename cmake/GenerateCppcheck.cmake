# - Generate a cppcheck documentation for a project.
# The function GENERATE_CPPCHECK is provided to create a "cppcheck" target that
# performs static code analysis using the cppcheck utility program.
#
# GENERATE_CPPCHECK(SOURCES <sources to check...>
#                   [SUPPRESSION_FILE <file>]
#                   [ENABLE_IDS <id...>]
#                   [TARGET_NAME <name>]
#                   [INCLUDES <dir...>]
#                   [INLINE_SUPPRESSION]
#                   [SYSTEM_INCLUDE_WARNING_SUPPRESSION]
#                   [HTML_REPORT])
#
# @brief Generates a target "cppcheck" that executes cppcheck on the specified 
# sources.
#
# @details 
# Cppcheck will be executed with CMAKE_CURRENT_SOURCE_DIR as working directory.
# This function can always be called, even if no cppcheck was found. Then no
# target is created.
#
# @param[in]    SOURCES
#
# @param[in]    SUPPRESSION_FILE   (required)  may be give additionally to 
# specify suppressions for cppcheck. The sources mentioned in the suppression 
# file must be in the same format like given for SOURCES. This means if you 
# specified them relative to CMAKE_CURRENT_SOURCE_DIR, then the same relative 
# paths must be used in the suppression file.
#
# @param[in]    ENABLE_IDS  allows to specify which additional cppcheck 
# check ids to execute, e.g. all or style. They are combined with AND.
#
# @param[in]    TARGET_NAME a different name for the generated check target can 
# be specified. This is useful if several calles to this function are made in one
# CMake project, as otherwise the target names collide.
#
# @param[in]    INCLUDES    additional include directories for the cppcheck 
# program can be given.
#
# @param[in]    INLINE_SUPPRESSION  is set, cppcheck inline-suppression comments 
# are parsed.
#
# @param[in]    SYSTEM_INCLUDE_WARNING_SUPPRESSION
#
# @param[in]    HTML_REPORT
#
# Copyright (C) 2011 by Johannes Wienke <jwienke at techfak dot uni-bielefeld dot de>
# Copyright (C) 2016 by Lionel Untereiner <untereiner dot lionel at gmail dot com>
#
# This program is free software; you can redistribute it
# and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation;
# either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

GET_FILENAME_COMPONENT(GENERATE_CPPCHECK_MODULE_DIR ${CMAKE_CURRENT_LIST_FILE} PATH)

FIND_PACKAGE(Cppcheck)

FUNCTION(GENERATE_CPPCHECK)

    IF(CPPCHECK_FOUND)

        set(options INLINE_SUPPRESSION SYSTEM_INCLUDE_WARNING_SUPPRESSION HTML_REPORT)
        set(oneValueArgs SOURCES SUPPRESSION_FILE ENABLE_IDS TARGET_NAME INCLUDES PLATFORM_FLAGS)
        set(multiValueArgs "")
    
        cmake_parse_arguments(ARG "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
        
        SET(TARGET_NAME "cppcheck")
        SET(TARGET_NAME_SUFFIX "")
        # parse target name
        LIST(LENGTH ARG_TARGET_NAME TARGET_NAME_LENGTH)
        IF(${TARGET_NAME_LENGTH} EQUAL 1)
            SET(TARGET_NAME ${ARG_TARGET_NAME})
            SET(TARGET_NAME_SUFFIX "-${ARG_TARGET_NAME}")
        ENDIF()
        
        SET(CPPCHECK_CHECKFILE "${CMAKE_BINARY_DIR}/cppcheck-files${TARGET_NAME_SUFFIX}")
        SET(CPPCHECK_REPORT_FILE "${CMAKE_BINARY_DIR}/cppcheck-report${TARGET_NAME_SUFFIX}.xml")
        SET(CPPCHECK_WRAPPER_SCRIPT "${CMAKE_BINARY_DIR}/cppcheck${TARGET_NAME_SUFFIX}.cmake")
    
        # write a list file containing all sources to check for the call to
        # cppcheck
        FILE(WRITE "${CPPCHECK_CHECKFILE}" "")

        MESSAGE(STATUS "Find in ${ARG_SOURCES}\n")
        FILE( GLOB_RECURSE SOURCE_FILES  ${ARG_SOURCES} "${ARG_SOURCES}/*.[ch]" "${ARG_SOURCES}/*.[ch]pp")

        FOREACH(SOURCE ${SOURCE_FILES})
            # MESSAGE(STATUS "Using sources file ${SOURCE}")
            FILE(APPEND "${CPPCHECK_CHECKFILE}" "${SOURCE}\n")
        ENDFOREACH()
        
        # prepare a cmake wrapper to write the stderr output of cppcheck to
        # the result file
        
        # suppression argument
        LIST(LENGTH ARG_SUPPRESSION_FILE SUPPRESSION_FILE_LENGTH)
        IF(${SUPPRESSION_FILE_LENGTH} EQUAL 1)
            GET_FILENAME_COMPONENT(ABS "${ARG_SUPPRESSION_FILE}" ABSOLUTE)
            MESSAGE(STATUS "Using suppression file ${ABS}")
            SET(SUPPRESSION_ARGUMENT --suppressions)
            SET(SUPPRESSION_FILE "\"${ABS}\"")
        ENDIF()

        # inline suppressions
        SET(INLINE_ARG)
        IF(ARG_INLINE_SUPPRESSION)
            SET(INLINE_ARG "--inline-suppr")
        ENDIF()

        SET(SYSTEM_INCLUDE_WARNING_SUPPRESSION_ARG)
        IF(ARG_SYSTEM_INCLUDE_WARNING_SUPPRESSION)
            SET(SYSTEM_INCLUDE_WARNING_SUPPRESSION_ARG "--suppress=missingIncludeSystem")
        ENDIF()

        # includes
        SET(INCLUDE_ARGUMENTS "")
        FOREACH(INCLUDE ${ARG_INCLUDES})
            SET(INCLUDE_ARGUMENTS "${INCLUDE_ARGUMENTS} \"-I${INCLUDE}\"")
        ENDFOREACH()
        
        # enabled ids
        SET(ID_LIST "")
        FOREACH(ID ${ARG_ENABLE_IDS})
            SET(ID_LIST "${ID_LIST},${ID}")
        ENDFOREACH()
        IF(ID_LIST)
            STRING(LENGTH ${ID_LIST} LIST_LENGTH)
            MATH(EXPR FINAL_LIST_LENGTH "${LIST_LENGTH} - 1")
            STRING(SUBSTRING ${ID_LIST} 1 ${FINAL_LIST_LENGTH} FINAL_ID_LIST)
            SET(IDS_ARGUMENT "\"--enable=${FINAL_ID_LIST}\"")
        ELSE()
            SET(IDS_ARGUMENT "")
        ENDIF()

        FILE(WRITE ${CPPCHECK_WRAPPER_SCRIPT}
"
EXECUTE_PROCESS(COMMAND \"${CPPCHECK_EXECUTABLE}\" ${INCLUDE_ARGUMENTS} ${SUPPRESSION_ARGUMENT} ${SUPPRESSION_FILE} ${INLINE_ARG} ${IDS_ARGUMENT} ${SYSTEM_INCLUDE_WARNING_SUPPRESSION_ARG} --inconclusive --force --xml-version=2  --language=c++ ${ARG_PLATFORM_FLAGS} \"--file-list=${CPPCHECK_CHECKFILE}\"
                RESULT_VARIABLE CPPCHECK_EXIT_CODE
                ERROR_VARIABLE ERROR_OUT
                WORKING_DIRECTORY \"${CMAKE_CURRENT_SOURCE_DIR}\")
IF(NOT CPPCHECK_EXIT_CODE EQUAL 0)
    MESSAGE(FATAL_ERROR \"Error executing cppcheck for target ${TARGET}, return code: \${CPPCHECK_EXIT_CODE}\")
ENDIF()
# IF(ERROR_OUT)
#     MESSAGE(\"Detected errors:\\n\${ERROR_OUT}\")
# ENDIF()
FILE(WRITE \"${CPPCHECK_REPORT_FILE}\" \"\${ERROR_OUT}\")
"
            )

        IF(ARG_HTML_REPORT)
            SET(CPPCHECK_REPORT_DIR "${CMAKE_BINARY_DIR}/cppcheck-report${TARGET_NAME_SUFFIX}")
            FILE(APPEND ${CPPCHECK_WRAPPER_SCRIPT}
"
EXECUTE_PROCESS(COMMAND \"${CPPCHECK_EXECUTABLE}-htmlreport\" \"--file=${CPPCHECK_REPORT_FILE}\" \"--title=CGoGN\" \"--source-dir=${ARG_SOURCES}\" \"--report-dir=${CPPCHECK_REPORT_DIR}\")
"
                )
        ENDIF()
            
        ADD_CUSTOM_TARGET(${TARGET_NAME} ${CMAKE_COMMAND} -P "${CPPCHECK_WRAPPER_SCRIPT}"
                          COMMENT "Generating cppcheck result ${TARGET_NAME}")

    ENDIF()

ENDFUNCTION()
 