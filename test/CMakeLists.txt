# Want test modules in the build directory.
#unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(gather_test gather_test.F90)
target_link_libraries(gather_test common MPI::MPI_Fortran)
add_test(gather ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./gather_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(gather PROPERTIES PROCESSORS 4)

add_executable(scatter_test scatter_test.F90)
target_link_libraries(scatter_test common MPI::MPI_Fortran)
add_test(scatter ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./scatter_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(scatter PROPERTIES PROCESSORS 4)

add_executable(localize_test localize_test.F90)
target_link_libraries(localize_test common MPI::MPI_Fortran)
add_test(localize ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./localize_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(localize PROPERTIES PROCESSORS 4)

add_executable(collate_test collate_test.F90)
target_link_libraries(collate_test common MPI::MPI_Fortran)
add_test(collate ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./collate_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(collate PROPERTIES PROCESSORS 4)

add_executable(distribute_test distribute_test.F90)
target_link_libraries(distribute_test common MPI::MPI_Fortran)
add_test(distribute ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./distribute_test ${MPIEXEC_POSTFLAGS})
set_tests_properties(distribute PROPERTIES PROCESSORS 4)
