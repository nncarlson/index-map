cmake_minimum_required(VERSION 3.19)
#set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

project(map_index Fortran)

option(USE_CAF "Build a version based on Fortran coarrays instead of MPI" OFF)

if(USE_CAF)
  add_compile_definitions(USE_CAF)
  add_compile_options($<$<Fortran_COMPILER_ID:NAG,Intel>:-coarray>)
  add_link_options($<$<Fortran_COMPILER_ID:NAG,Intel>:-coarray>)
else()
  find_package(MPI REQUIRED)
endif()
add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>)

#add_compile_definitions($<$<COMPILE_LANGUAGE:Fortran>:NDEBUG>)

include(CTest)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(example)
