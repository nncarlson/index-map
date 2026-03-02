cmake_minimum_required(VERSION 3.19)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

project(map_index Fortran)

option(USE_CAF "Build a version based on Fortran coarrays instead of MPI" OFF)
option(BUILD_HTML "Build Sphinx HTML documentation" OFF)

find_package(Python REQUIRED COMPONENTS Interpreter)
find_program(FYPP_EXECUTABLE fypp)
if(NOT FYPP_EXECUTABLE)
  message(FATAL_ERROR
      "Could not find the 'fypp' preprocessor.\n"
      "Install it with 'python3 -m pip install fypp' and rerun CMake.")
endif()
include(Fypp)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build" FORCE)
endif()

add_compile_definitions("$<$<NOT:$<CONFIG:Debug>>:NDEBUG>")

if(USE_CAF)
  add_compile_definitions(USE_CAF)
  add_compile_options($<$<Fortran_COMPILER_ID:NAG,Intel,IntelLLVM>:-coarray>)
  add_link_options($<$<Fortran_COMPILER_ID:NAG,Intel,IntelLLVM>:-coarray>)
else()
  find_package(MPI REQUIRED)
endif()

add_compile_options($<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffree-line-length-none>)

# Set missing CMake defaults for the NAG Fortran compiler
if(CMAKE_Fortran_COMPILER_ID STREQUAL NAG)
  set(CMAKE_Fortran_FLAGS_DEBUG "-u -O0 -gline -C -nan")
  set(CMAKE_Fortran_FLAGS_RELEASE "-u -O3")
  set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELEASE} -gline")
endif()

include(CTest)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(example)

if(BUILD_HTML)
  find_program(SPHINX_EXECUTABLE NAMES sphinx-build
    DOC "Path to sphinx-build executable"
  )
  if(NOT SPHINX_EXECUTABLE)
    message(FATAL_ERROR
        "Could not find 'sphinx-build'.\n"
        "Install Sphinx and rerun CMake.")
  endif()
  add_subdirectory(doc/sphinx)
endif()
