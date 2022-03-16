cmake_minimum_required(VERSION 3.19)

project(map_index Fortran)

find_package(MPI REQUIRED)

add_compile_definitions($<$<COMPILE_LANGUAGE:Fortran>:NDEBUG>)

include(CTest)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(example)
