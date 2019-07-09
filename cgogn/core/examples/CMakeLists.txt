cmake_minimum_required(VERSION 3.7.2 FATAL_ERROR)

project(cgogn_core_examples
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_io REQUIRED)

add_executable(core_test core_test.cpp)
target_link_libraries(core_test cgogn::io cgogn::core)

set_target_properties(core_test PROPERTIES FOLDER examples/core)
