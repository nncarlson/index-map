# Performance Results

## Hash 93b316

Update: spot checks with ec866b show same timings.

* index-map source hash 93b316 with some modifications (noted)
* All MPICH versions built with `--with-device=ch3`
* Timings (usec) for `disk-fv-mpi` modified with `NZ=257` using 8 processes
  on 12-core AMD Threadripper 2920X. Test run 7 times; 2 shortest and 2 longest
  times discarded and the remaining 3 times averaged.
* Two code versions considered:
  - Version 1: Compact with vector subscript (original)
    ```fortran
    u_local(j) = u_prev(j) + c*(sum(u_prev(cnhbr_local(:,j))) - 4*u_prev(j))
    ```
  - Version 2: Explicit do-loop (modified)
    ```fortran
    s = 4*u_prev(j)
    do k = 1, size(cnhbr_local,1)
      s = s - u_prev(cnhbr_local(k,j))
    end do
    u_local(j) = u_prev(j) - c*s
    ```
#### General Observations

* MPICH seems to perform slightly better than OpenMPI
* MPICH performance seems to have slightly degraded with each subsequent version
* Intel and GFortran (with higher optimization) are agnostic with respect to
  code version.
* NAG requires the explicit do-loop version with higher optimization to
  achieve comparable results. If not, performance is very poor.
* Best results: GFortran (~17), NAG (~20), Intel (~23)
### Intel Compiler

intel-19.1.0    |V1| -O2  | -O3  |V2| -O2  | -O3  |
----------------|--|------|------|--|------|------|
mpich-3.3.2     |  | 20.3 | 23.4 |  | 23.4 | 22.0 |
mpich-3.4.3     |  | 20.8 | 23.0 |  | 23.9 | 22.1 |
mpich-4.0.1     |  | 20.7 | 20.9 |  | 24.0 | 23.7 |
openmpi-4.1.2   |  | 21.3 | 21.6 |  | 24.8 | 22.2 |

* Best results are with the original compact code using default `-O2`
  optimization.
* Higher optimization `-O3` provides small benefit for the explicit
  do-loop version, but performs worse for the compact code.

intel-2021.5.0  |V1| -O2  | -O3  |V2| -O2  | -O3  |
----------------|--|------|------|--|------|------|
mpich-3.3.2     |  | 23.4 | 23.4 |  | 23.3 | 23.4 |
mpich-3.4.3     |  | 23.4 | 23.5 |  | 23.6 | 23.4 |
mpich-4.0.1     |  | 23.6 | 23.8 |  | 23.6 | 23.8 |
openmpi-4.1.2   |  | 24.8 | 24.6 |  | 21.5 | 24.5 |

* Higher optimization provides no benefit whatsoever.
* No difference between the two code versions.
* Noticeably worse performance than with Intel 19.1

### GFortran Compiler

gfortran-11.2.0 |V1| -O2  | -O3  |V2| -O2  | -O3  |
----------------|--|------|------|--|------|------|
mpich-3.3.2     |  | 20.3 | 17.1 |  | 17.0 | 17.0 |
mpich-3.4.3     |  | 20.7 | 17.5 |  | 17.4 | 17.6 |
mpich-4.0.1     |  | 20.9 | 17.7 |  | 17.8 | 17.9 |
openmpi-4.1.2   |  | 17.9 | 18.1 |  | 21.2 | 18.0 |

* Significantly better performance than with the other compilers.
* Compact code performs as well explicit do-loop code but only
  when higher optimization is used.
* Why performance improves under higher optimization for the
  explicit do-loop version with openmpi but not mpich is a
  mystery.

### NAG Compiler

nagfor-7.0.7026 |V1| -O2  | -O3 | -O4 |V2| -O2  | -O3  | -O4  |
----------------|--|------|-----|-----|--|------|------|------|
mpich-3.3.2     |  | 92.0 | 114 | 114 |  | 27.9 | 19.8 | 20.0 |
mpich-3.4.3     |  | 92.0 | 114 | 114 |  | 28.0 | 20.2 | 20.2 |
mpich-4.0.1     |  | 93.6 | 114 | 114 |  | 28.1 | 20.6 | 20.4 |
openmpi-4.1.2   |  | 93.3 | 115 | 115 |  | 28.7 | 20.9 | 21.2 |

mpich-4.0.1     |V1| -O2  | -O3 | -O4 |V2| -O2  | -O3  | -O4  |
----------------|--|------|-----|-----|--|------|------|------|
nagfor-7.0.7026 |  | 93.6 | 114 | 114 |  | 28.1 | 20.6 | 20.4 |
nagfor-7.1.7105 |  | 91.2 | 112 | 112 |  | 28.1 | 20.5 | 20.3 |

openmpi-4.1.2   |V1| -O2  | -O3 | -O4 |V2| -O2  | -O3  | -O4  |
----------------|--|------|-----|-----|--|------|------|------|
nagfor-7.0.7026 |  | 93.3 | 115 | 115 |  | 28.7 | 20.9 | 21.2 |
nagfor-7.1.7105 |  | 91.5 | 113 | 113 |  | 28.6 | 20.9 | 21.2 |

* Performance of the compact code is very poor (factor 3.5-4) under any
  optimization level, and actually gets significantly worse for the higher
  levels.
* With the explicit do-loop code, a higher `-O3` optimization gives very
  much improved results, but no further gain with an even higher level.
* Results from version 7.0 and 7.1 are pretty much the same.
