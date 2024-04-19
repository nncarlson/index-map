### Timings for the `disk-fv` Example (18 Apr 2024, 441e52b)

* [Test configuration](./test-configuration.md)
* 12-core AMD Threadripper 2920X cpu
* Values are average µsec per time step
* DNF = program segfaulted or was terminated after a very long time.

#### Using `NZ=513`

| parallel  |   1 |    2 |    4 |    6 |    8 |   12 |
|-----------|----:|-----:|-----:|-----:|-----:|-----:|
| NAG-MPI   | 547 |  281 |  140 |   95 |   72 |   52 |
| NAG-CAF   | 341 |  157 |  102 |   82 |   64 |   59 |
| GNU-MPI   | 322 |  167 |   83 |   69 |   61 |   45 |
| GNU-CAF   | 455 | 1245 | 3834 | 3841 | 3755 | 4062 |
| ifort-MPI | 678 |  385 |  223 |  150 |  112 |   78 |
| ifx-MPI   | 693 |  384 |  215 |  144 |  110 |   77 |

#### Using the default `NZ=257`

|        | NAG | GNU | ifort | ifx |
|--------|-----|-----|-------|-----|
| serial | 131 | 100 |   113 | 167 |

| parallel  |   1 |   2 |    4 |    6 |    8 |   12 |
|-----------|----:|----:|-----:|-----:|-----:|-----:|
| NAG-MPI   | 131 |  67 |   36 |   25 |   20 |   15 |
| NAG-CAF   |  74 |  46 |   33 |   31 |   30 |   31 |
| GNU-MPI   |  74 |  39 |   22 |   15 |   12 |   10 |
| GNU-CAF   | 109 | 596 | 1975 | 2006 | 1934 | 2123 |
| ifort-MPI | 170 |  98 |   50 |   32 |   25 |   18 |
| ifort-CAF | 159 | DNF |  DNF |  DNF |  DNF |  DNF |
| ifx-MPI   | 172 |  96 |   46 |   31 |   25 |   18 |
| ifx-CAF   | 160 | DNF |  DNF |  DNF |  DNF |  DNF |
| ifx-MPICH | 175 |  92 |   48 |   30 |   23 |   17 |

#### Using `NZ=127`

| parallel  |   1 |   2 |    4 |    8 |
|-----------|----:|----:|-----:|-----:|
| NAG-MPI   |  31 |    17 |    10 |     7 |
| NAG-CAF   |  18 |    16 |    19 |    22 |
| GNU-MPI   |  18 |    10 |     7 |     5 |
| GNU-CAF   |  27 |   376 |  1097 |  1127 |
| ifort-MPI |  39 |    22 |    13 |     9 |
| ifort-CAF |  40 |  1257 |  2994 |  3352 |
| ifx-MPI   |  38 |    21 |    12 |     8 |
| ifx-CAF   |  39 |  1230 |  2952 |  3290 |
| ifx-MPICH |  38 |    20 |    12 |     8 |
