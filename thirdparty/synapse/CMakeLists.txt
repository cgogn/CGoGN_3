project(synapse LANGUAGES C CXX)

set(HEADER_FILES
	include/boost/synapse/block.hpp
	include/boost/synapse/blocker.hpp
	include/boost/synapse/connect.hpp
	include/boost/synapse/connection.hpp
	include/boost/synapse/emit.hpp
	include/boost/synapse/properties.hpp
	include/boost/synapse/signal_traits.hpp
	include/boost/synapse/thread_local_queue.hpp
	include/boost/synapse/translate.hpp
	include/boost/synapse/dep/assert.hpp
	include/boost/synapse/dep/functional.hpp
	include/boost/synapse/dep/smart_ptr.hpp
	include/boost/synapse/dep/thread_local.hpp
	include/boost/synapse/synapse_detail/common.hpp
	include/boost/synapse/synapse_detail/config.hpp
	include/boost/synapse/synapse_detail/weak_store.hpp
)

set(SOURCE_FILES
	src/block.cpp
	src/connect.cpp
	src/thread_local_queue.cpp
)

add_definitions(-DBOOST_SYNAPSE_USE_STD_SMART_PTR -DBOOST_SYNAPSE_USE_STD_FUNCTION -DBOOST_SYNAPSE_ASSERT)

add_library(${PROJECT_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES})
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")

target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>/include
	$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cgogn/thirdparty>
)

install(FILES ${HEADER_FILES} DESTINATION "include/cgogn/thirdparty")

cgogn_create_package("${CMAKE_CURRENT_SOURCE_DIR}")

set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER thirdparty)