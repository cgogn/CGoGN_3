cmake_minimum_required(VERSION 3.7.2 FATAL_ERROR)

project(cgogn_io_examples
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)
find_package(cgogn_io REQUIRED)

add_executable(io_obj_test io_obj_test.cpp)
target_link_libraries(io_obj_test cgogn::io cgogn::core)

set_target_properties(io_obj_test PROPERTIES FOLDER examples/io)
