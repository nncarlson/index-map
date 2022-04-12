### Timings for the `disk-fv` Example (11 Apr 2022, 34d22ed)

* [Test configuration](./test-configuration.md)
* 12-core AMD Threadripper 2920X cpu
* Values are average µsec per time step
* DNF = program terminated or segfaulted after a very long time.

#### Using the default `NZ=257`

|        | NAG | GNU | Intel |
|--------|-----|-----|-------|
| serial | 131 |  75 |   102 |

| parallel  |   1 |   2 |    4 |    6 |    8 |   12 |
|-----------|----:|----:|-----:|-----:|-----:|-----:|
| NAG-MPI   | 137 |  70 |   37 |   26 |   20 |   16 |
| NAG-CAF   |  87 |  51 |   34 |   29 |   28 |   29 |
| GNU-MPI   |  76 |  39 |   22 |   15 |   13 |   10 |
| GNU-CAF   | 110 | 434 | 1538 | 1485 | 1454 | 1608 |
| Intel-MPI | 134 |  70 |   37 |   27 |   21 |   17 |
| Intel-CAF | 171 | DNF |  DNF |  DNF |  DNF |  DNF |

#### Using `NZ=127`

| parallel  |   1 |   2 |    4 |    8 |
|-----------|----:|----:|-----:|-----:|
| NAG-MPI   |  32 |    20 |    15 |    11 |
| NAG-CAF   |  21 |    20 |    20 |    21 |
| GNU-MPI   |  19 |    14 |    11 |     8 |
| GNU-CAF   |  26 |   244 |   815 |   820 |
| Intel-MPI |  50 |    33 |    24 |    19 |
| Intel-CAF |  44 | 25000 | 43000 | 79000 |
