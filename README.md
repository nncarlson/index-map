## Index Map

The `index_map_type` Fortran module provides core capabilities that support
use of distributed arrays in MPI-based SPMD programs through the `index_map`
derived type that describes the mapping of an array's index set to processes.
The mapping allows for overlap between processes, and provides collective
array gather procedures and scatter reduction procedures associated with that
overlap.

Documentation for the module is in the `doc` directory.

Usage examples are in the `example` directory.
