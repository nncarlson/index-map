set(SRC
  f90_assert.F90
  integer_set_type.F90
  integer_map_type.F90
  index_map_type.F90
  index_map_type-gather_offp_impl.F90
  index_map_type-scatter_offp_impl.F90
  index_map_type-distribute_impl.F90
  index_map_type-collate_impl.F90
  index_map_type-localize_impl.F90
)

add_library(common OBJECT ${SRC})
target_link_libraries(common MPI::MPI_Fortran)
target_include_directories(common PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)

# We are using MPI's Fortran "mpi" module and call the MPI functions
# with different data types. Some compilers will detect this and report
# compilation errors. Here we add compilation flags that tell the compiler
# to relax and not worry about the apparent contradictory usages.
# NB: This issue goes away if we were to use the "mpi_f08" module instead.

set(MPI_USERS
  index_map_type.F90
  index_map_type-gather_offp_impl.F90
  index_map_type-scatter_offp_impl.F90
  index_map_type-distribute_impl.F90
  index_map_type-collate_impl.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES NAG)
  set_source_files_properties(${MPI_USERS}
      PROPERTIES COMPILE_OPTIONS "-mismatch;-w")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  set_source_files_properties(${MPI_USERS}
      PROPERTIES COMPILE_OPTIONS "-fallow-argument-mismatch;-w")
endif()
