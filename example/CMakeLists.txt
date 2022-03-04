# Want test modules in the build directory.
#unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(disk-fv-serial disk-fv-serial.f90)

add_executable(disk-fv-mpi disk-fv-mpi.f90)
target_link_libraries(disk-fv-mpi common MPI::MPI_Fortran)
