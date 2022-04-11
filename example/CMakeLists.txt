# Want test modules in the build directory.
#unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(disk-fv-serial disk-fv-serial.F90)

add_executable(disk-fv-parallel disk-fv-parallel.F90)
target_link_libraries(disk-fv-parallel common)

add_executable(disk-fem-serial disk-fem-serial.F90)

add_executable(disk-fem-parallel disk-fem-parallel.F90)
target_link_libraries(disk-fem-parallel common)
