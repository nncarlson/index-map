## Index Map

The `index_map_type` Fortran module provides core capabilities that support
the use of distributed arrays in SPMD parallel programs through the `index_map`
derived type, which describes the mapping of an array's index set to parallel
processes. The mapping allows for overlap between processes, and provides
collective array gather procedures and scatter-reduction procedures associated
with that overlap.

There are two implementations of the module: one using MPI located in
`src/mpi`, and one using Fortran coarrays located in `src/caf`. The
interfaces are virtually the same, making it relatively simple for client
code to switch between the two versions.

* The MPI version is stable and in production use
(see [Truchas](https://gitlab.com/truchas/truchas)).

* The coarray version is beta quality. It is fully tested, but its current
performance is highly variable. Performance of the example programs is
comparable to the MPI version when using the NAG Fortran compiler and its
shared memory-only implementation of coarrays. However performance is very
poor when using GFortran / OpenCoarrays, and worse still with the Intel
Fortran compiler and its implementation of coarrays. Links to performance
results can be found in this [README](./doc/README.md). This is a work in
progress.

### Documentation

Documentation for the module is in the file
[`doc/index_map_type.md`](./doc/index_map_type.md).

The `doc` directory [README](./doc/README.md) has links to some additional
information.

### Examples
The `example` directory contains several example programs that illustrate
the usage of the `index_map_type` module. These programs are automatically
compiled as part of the build process. See its [README](./example/README.md)
file for details on the examples and how to run them.

### Compiling

First clone the repository:
```
$ git clone https://github.com/nncarlson/index-map.git
$ cd index-map
```
The project uses CMake (version 3.19 or later) to compile the project.
You may need to set your `FC` environment variable to the path of your
Fortran compiler before running `cmake` to ensure that CMake finds the
correct compiler.

#### MPI version
```
$ mkdir build
$ cd build
$ cmake .. -DCMAKE_Fortran_FLAGS="-O3 -DNDEBUG"
$ make
```
Set your desired compiler flags using `CMAKE_Fortran_FLAGS`. Define the
preprocessor macro `NDEBUG` (as shown) to disable runtime assertion checks.

If CMake does not automatically find the MPI installation that matches
your Fortran compiler, you can help it by setting the environment variable
`MPI_HOME` to the root directory of the MPI installation. Note that the MPI
Fortran interface (the `mpi` module) must have been compiled with the same
Fortran compiler being used to compile this project.

#### Fortran coarray version
This CMake setup understands how to build the coarray version when using one
of these Fortran compilers:
* NAG 7.1 or later with its built-in coarray support.
* Intel with its built-in coarray support. The companion Intel MPI package
  must be installed and Intel's setup script run to configure your environment.
* GFortran with [OpenCoarrays](https://github.com/sourceryinstitute/opencoarrays).
  OpenCoarrays supplies the implementation of coarrays used by the gfortran
  compiler. Be sure the `bin` directory of the opencoarrays installation is in
  your path so that the compiler wrapper `caf` and runner `cafrun` can be found.
  Set `FC=caf` before running cmake.

```
$ mkdir build
$ cd build
$ cmake .. -DUSE_CAF=Yes -DCMAKE_Fortran_FLAGS="-O3 -DNDEBUG"
$ make
```
### Installation
The CMake setup does not currently install a library. It will merely build
the unit test programs and example programs. The current thinking is that
a user will simply incorporate the source code for the `index_map_type`
module directly into the source for their application. This could change in
the future.

### Testing
From the `build` directory simply run `ctest`:
```
$ ctest
Test project /home/user/index-map/build
    Start 1: localize
1/5 Test #1: localize .........................   Passed    0.07 sec
    Start 2: gather
2/5 Test #2: gather ...........................   Passed    0.07 sec
    Start 3: scatter
3/5 Test #3: scatter ..........................   Passed    0.08 sec
    Start 4: collate
4/5 Test #4: collate ..........................   Passed    0.07 sec
    Start 5: distribute
5/5 Test #5: distribute .......................   Passed    0.08 sec

100% tests passed, 0 tests failed out of 5

Total Test time (real) =   0.37 sec
```
Each test is actually a collection of sub-tests (about 100 in all);
run `ctest -V` to get the verbose output. All tests should pass.

Note: These parallel tests use 4 processes. Without additional options /
environment settings, MPI may refuse to run "over-subscribed" when testing
on a system with fewer than 4 cores or hardware threads.  

### License
This project is distributed under the terms of the MIT license.
See [LICENSE.md](./LICENSE.md) for details.
