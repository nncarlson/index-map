### Timings for the `disk-fem` Example (18 Apr 2024, 441e52b)

* [Test configuration](./test-configuration.md)
* 12-core AMD Threadripper 2920X cpu
* Values are average µsec per time step
* DNF = program segfaulted or was terminated after a very long time.

#### Using the default `NZ=257`

|        | NAG | GNU | ifort |
|--------|-----|-----|-------|
| serial | 288 | 178 |   200 |

| | 1 | 2 | 4 | 6 | 8 | 12 |
|-|---|---|---|---|---|----|
| NAG-MPI   | 268 |  140 |   73 |   50 |   40 |   29 |
| NAG-CAF   | 221 |  130 |   83 |   71 |   67 |   68 |
| GNU-MPI   | 179 |   93 |   50 |   35 |   28 |   21 |
| GNU-CAF   | 184 | 1687 | 4550 | 4880 | 5412 | 5793 |
| ifort-MPI | 200 |  105 |   55 |   38 |   31 |   24 |
| ifort-CAF | 200 |  DNF |  DNF |  DNF |  DNF |  DNF |

#### Using `NZ=101`

| | 1 | 2 | 4 | 8 |
|-|---|---|---|---|
| ifort-MPI | 31 |   19 |   12 |    9 |
| ifort-CAF | 31 | 1372 | 2616 | 3003 |
