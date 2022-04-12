## Example Programs

This directory contains example programs that illustrate the use of the
`index_map_type` Fortran module.

The compiled executables require no input. The serial executables are run like
any normal executable. Running the parallel executables varies. For example:
* MPI version:
  ```
  $ mpirun -np <N> ./disk-fv-parallel
  ```
* Coarray version:
  - NAG:
    ```
    $ NAGFORTRAN_NUM_IMAGES=<N> ./disk-fv-parallel
    ```
  - Intel:
    ```
    $ FOR_COARRAY_NUM_IMAGES=<N> ./disk-fv-parallel
    ```
  - GFortran:
    ```
    $ cafrun -n <N> ./disk-fv-parallel
    ```

### Finite Volume Solution of the Heat Equation on a Unit Disk

The program `disk-fv-parallel.F90` solves the heat equation $`u_t = \Delta u`$
on the unit disk subject to zero boundary conditions. It uses a simple finite
volume discretization on a regular 2D Cartesian mesh that covers the unit
disk. Only those cells whose center is contained in the unit disk are
included in the problem. The unknowns are cell-centered values of $`u`$.
Simple first-order forward Euler time stepping is used to solve from a
uniform initial condition to a final time. Domain decomposition is used to
parallelize the computation. The cells are numbered according to the usual
order of the elements of a rank-2 array (but only for those included in
the problem) and then partitioned into approximately equal blocks, one per
process (i.e., MPI rank or coarray image). There is a single indirect indexing
array that maps a cell index to the cell indices of its neighboring cells. For
comparison, the program `disk-fv-serial.F90` is a serial version of the solver.
Aside from the initial setup of the parallel decomposition of the problem, the
only difference between the parallel solver and the serial solver is a single
call to `index_map%gather_offp` to update the off-process elements of the local
unknown array with values from their corresponding on-process elements on
neighboring processes -- a parallel halo exchange operation.

### Finite Element Solution of the Heat Equation on a Unit Disk

The program `disk-fem-parallel.F90` solves the heat equation $`u_t = \Delta u`$
on the unit disk subject to zero boundary conditions. It uses a Galerkin finite
element discretization with bilinear elements on a regular 2D Cartesian mesh
that covers the unit disk. Only those cells with at least one node contained
in the interior of the unit disk are included in the problem. The unknowns are
nodal values of $`u`$ for nodes contained in the interior of the unit disk.
In order to avoid complexities associated with the need to solve a linear
system, a lumped mass matrix and simple first-order forward Euler time stepping
are used to solve from a uniform initial condition to a final time. Domain
decomposition is used to parallelize the computation. The cells are numbered
according to the usual order of the elements of a rank-2 array (but only for
those included in the problem) and then partitioned into approximately equal
blocks, one per process (i.e., MPI rank or coarray image). The nodes are
similarly numbered, and are partitioned compatibly with the cells: a node is
assigned to a partition of one of its adjacent cells. There is a single
indirect indexing array that maps a cell index to the indices of its adjacent
nodes. For comparison, the program `disk-fem-serial.F90` is a serial version
of the solver. Aside from the initial setup of the parallel decomposition of
the problem, the only difference in the parallel solver are two parallel halo
exchange calls: a call to `index_map%gather_offp` at the beginning of each
time step to update the off-process elements of the local unknown array, and
a later call to `index_map%scatter_offp_sum` to complete the FE assembly of
the right-hand-side Laplacian term.
