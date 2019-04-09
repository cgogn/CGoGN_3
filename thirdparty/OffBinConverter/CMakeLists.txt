project(off_ascii2bin
	LANGUAGES CXX
)

find_package(cgogn_core REQUIRED)

add_executable(off_ascii2bin off_ascii2bin.cpp)
target_link_libraries(off_ascii2bin cgogn::core)

if(NOT MSVC)
	if (CGOGN_CPP_STD STREQUAL "11")
		target_compile_options(off_ascii2bin PUBLIC "-std=c++11")
	elseif(CGOGN_CPP_STD STREQUAL "14")
		target_compile_options(off_ascii2bin PUBLIC "-std=c++14")
	else(CGOGN_CPP_STD STREQUAL "11")
		target_compile_options(off_ascii2bin PUBLIC "-std=c++17")
	endif(CGOGN_CPP_STD STREQUAL "11")
endif(NOT MSVC)

set_target_properties(off_ascii2bin PROPERTIES FOLDER thirdparty)
