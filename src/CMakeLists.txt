set(SRC index_map_type.F90 index_map_type-gather_impl.F90)
add_library(common OBJECT ${SRC})
target_link_libraries(common MPI::MPI_Fortran)
target_include_directories(common PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
)
