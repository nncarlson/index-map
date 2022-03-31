## Index Map

The `index_map_type` Fortran module provides core capabilities that support
use of distributed arrays in SPMD parallel programs through the `index_map`
derived type that describes the mapping of an array's index set to processes.
The mapping allows for overlap between processes, and provides collective
array gather procedures and scatter reduction procedures associated with that
overlap.

There are two implementations of the module: one based on MPI located in
`src/mpi`, and one using Fortran coarrays located on `src/caf`. They have
essentially the same interfaces.

Note: The coarray version is alpha quality: it works, but its current
perfomance is very poor relative to the MPI version.

By default CMake will configure a build of the MPI version. Use the `cmake`
command line flag `-DUSE_CAF=yes` to build the coarray version instead.
Currently, the CMake setup only understands how to build the coarray version
when using the NAG compiler.

Documentation for the module is in the `doc` directory.

Usage examples are in the `example` directory.
