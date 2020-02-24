cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

project(cgogn_rendering_examples
        LANGUAGES CXX)

find_package(cgogn_core REQUIRED)
find_package(cgogn_ui REQUIRED)
find_package(cgogn_io REQUIRED)
find_package(cgogn_rendering REQUIRED)
find_package(cgogn_modeling REQUIRED)

set(APPS_LIST
#    modeling
#    tubular_mesh
#    deformation
#    mr_volume
    )

foreach(app ${APPS_LIST})
    apps_compilation(EXEC  "${app}"
        SRC  "${app}.cpp" 
        MODULES ui io rendering modeling)
endforeach()

