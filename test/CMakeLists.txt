# Want test modules in the build directory.
#unset(CMAKE_Fortran_MODULE_DIRECTORY)

add_executable(test-gather gather_test.F90)
target_link_libraries(test-gather common MPI::MPI_Fortran)

add_test(gather ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./test-gather ${CMAKE_CURRENT_SOURCE_DIR}/test-data/opencalc-B0-4 ${MPIEXEC_POSTFLAGS})
set_tests_properties(gather PROPERTIES PROCESSORS 4)
