## Some Coarray Performance Results (4 Apr 2022, cd26e73)

### A significant NAG performance issue
* Using the **serial** code `disk-fv-serial.F90` (NZ=257)
* Tests using nagfor build 7106 (same results with 7103)
* The `COMPACT_UPDATE` option uses this update statement
  ```fortran
  u(j) = u_prev(j) + c*(sum(u_prev(cnhbr(:,j))) - 4*u_prev(j))
  ```
  instead of the default explicit do loop:
  ```fortran
  tmp = 4*u_prev(j)
  do k = 1, size(cnhbr,1)
    tmp = tmp - u_prev(cnhbr(k,j))
  end do
  u(j) = u_prev(j) - c*tmp
  ```
* Run with `NAGFORTRAN_NUM_IMAGES=1`
* `-O4` gave same timings as `-O3`
* The code makes no use whatsoever of coarrays, so why adding the `-coarray` option
  makes any difference, let alone such a huge difference, is baffling. The timings
  are for just the computational loop; program startup and shutdown are not included.
* Without `-coarray`, it had already been observed that the compact code is far
  slower (4x) than the equivalent do loop version without `-coarray`; unlike
  gfortran and Intel, the NAG compiler is not able to optimize the compact code.
* But what is bizarre, is that adding `-coarray` (which one would think should
  have no impact) greatly reduces the time for the do loop version, and greatly
  increases the time for the compact code.

additional options        |  -O2 |  -O3 |
--------------------------|------|------|
(none)                    |  144 |  130 |
-DCOMPACT_UPDATE          |  647 |  855 |
-coarray                  |   88 |   73 |
-DCOMPACT_UPDATE -coarray | 1699 | 1704 |

### Parallel Tests

* Using the parallel code `disk-fv-parallel.F90` (NZ=257)
* Unless noted, the executables were **not** built using `-DCOMPACT_UPDATE`;
  that is the update step is implemented using a do loop. This is essential
  for decent NAG performance, and in other tests did not have a significant
  impact with either Intel or GFortran.

#### NAG

* nagfor build 7106, MPICH 4.0.1 (not used for CAF)
* `cmake -DCMAKE_Fortran_FLAGS="-O<2|3> -DNDEBUG" -DUSE_CAF=<yes|no>``
* MPI: run with `mpirun -np <n> -bind-to rr ...`
* CAF: run with `NAGFORTRAN_NUM_IMAGES=<n> ...`
* For some the time spent in computation is shown in (.); the balance of the
  time is spent in off-process gather communication.

|          | serial |   1 |    2    |    4    |    8    |
|----------|--------|-----|---------|---------|---------|
| MPI: -O2 |   145  | 175 | 89      | 48      | 26      |
| MPI: -O3 |   131  | 138 | 71 (70) | 37 (35) | 20 (18) |
| CAF: -O2 |   145  |  81 | 48 (43) | 33 (22) | 28 (12) |
| CAF: -O3 |   131  |  88 | 51      | 34      | 28      |

#### GFortran

* gfortran 11.2.0, mpich 4.0.1, opencoarrays 2.9.2-13-g235167d
* `cmake -DCMAKE_Fortran_FLAGS="-O<2|3> -DNDEBUG" -USE_CAF=<yes|no>`
* For CAF set `FC=caf` (instead of `gfortran`)
* MPI: run with `mpirun -np <n> -bind-to rr ...`
* CAF: run with `cafrun -n <n> -bind-to rr ...`
* For some the time spent in computation is noted in (.); the balance of the
  time is spent in the off-process communication.
* **These are atrocious times**

|          | serial |  1  |   2  |   4  |   8  |
|----------|--------|-----|------|------|------|
| MPI: -O2 |   114  | 109 |   56      |   30      |   17      |
| MPI: -O3 |    75  |  76 |   40 (38) |   23 (19) |   13 (10) |
| CAF: -O2 |   114  | 109 | 1260      | 4862      | 5029      |
| CAF: -O3 |    75  | 109 | 1258 (56) | 4850 (29) | 5025 (16) |


### Intel

Using ifort 2021.5.0

##### MPICH
* MPICH 4.0.1
* `mpirun -np <n> -bind-to rr ...`
 * Runtime segfault (after running for awhile) with CAF and mpich (4.0.1 and 3.3.2)

| MPICH | serial | 1 | 2 | 4 | 8 |
|----------|--------|---|---|---|---|
| -O2 |  97 | 146 | 74 | 39 | 21 |
| -O3 | 100 | 134 | 69 | 36 | 20 |

##### Intel MPI
*
      export I_MPI_FABRICS=shm:ofi
      export I_MPI_PIN_PROCESSOR_LIST=0-11
* `FOR_COARRAY_NUM_IMAGES=<n> disk-fv-parallel`
* Why is the serial time so much less than the 1-process MPI time?
  An extra call is made each time step than returns immediately;
  no MPI functions are called.
* Times are pretty much the same as MPICH, which is not surprising
  since Intel MPI is based on MPICH.

|          | serial | 1 | 2 | 4 | 8 |
|----------|--------|---|---|---|---|
| -O2 |  97 |  146 | 75 | 40 | 23 |
| -O3 | 100 |  134 | 69 | 37 | 21 |

##### Coarray

* Compile options: `-coarray -O<n>`
* Environment variables:
  ```
  export I_MPI_FABRICS=shm
  export I_MPI_DEVICE=shm
  export I_MPI_PIN_PROCESSOR_LIST=0-11
  ```
* 1-image executable segfaults during initialization due to a
  compiler/runtime [bug](https://community.intel.com/t5/Intel-Fortran-Compiler/Coarray-runtime-bug/m-p/1374152#M160881)
* Multi-image executables all segfault after a long time. Some tracebacks
  indicate it happens in the CAF runtime invoked from the off-process gather
  procedure.
* The following timings were obtained by reducing the size of the problem
  from NZ=257 to NZ=100.
* The time spent in the computation is in (.); the rest is the off-process
  communication call. **The results are atrocious**


|     | serial |  1  |      2      |      4      |      8      |
|-----|--------|-----|-------------|-------------|-------------|
| -O2 |    97  | --- | 14000 (7.6) | 30500 (4.7) | 58600 (2.3) |
| -O3 |   100  | --- | 14000 (11)  | 30100 (6)   | 57700 (3)   |
