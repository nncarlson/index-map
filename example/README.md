## Example Programs

This directory contains example programs that illustrate the use of the
`index_map_type` Fortran module.

### Finite Volume Solution of the Heat Equation on a Unit Disk

The program `disk-fv-mpi.f90` solves the heat equation $`u_t = \Delta u`$ on
the unit disk subject to zero boundary conditions. It uses a simple finite
volume discretization on a regular 2D Cartesian mesh that contains the unit
disk. Only those cells whose center is contained in the unit disk are
included in the problem. The unknowns are cell-centered values of $`u`$.
Simple first-order forward Euler time stepping is used to solve from a
uniform initial condition to a final time. Domain decomposition is used to
parallelize the computation. The cells are numbered according to the usual
order of the elements of a rank-2 array (but only for those included in
the problem) and then partitioned into approximately equal blocks, one per
MPI rank. There is a single indirect indexing array that maps a cell index
to the cell indices of its neighboring cells. For comparison, the program
`disk-fv-serial.F90` is a serial version of the solver. Aside from the initial
setup of the parallel decomposition of the problem, the only difference
between the MPI-parallel solver and the serial solver is a single call to
`index_map%gather_offp` to update the off-process elements of the local
unknown array with values from their corresponding on-process elements on
neighboring processes -- a parallel halo exchange operation.
