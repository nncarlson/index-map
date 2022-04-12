### Timings for the `disk-fem` Example (12 Apr 2022, 34d22ed)

* [Test configuration](./test-configuration.md)
* 12-core AMD Threadripper 2920X cpu
* Values are average µsec per time step
* DNF = program terminated or segfaulted after a very long time.

#### Using the default `NZ=257`

|        | NAG | GNU | Intel |
|--------|-----|-----|-------|
| serial | 322 | 177 |   200 |

| | 1 | 2 | 4 | 6 | 8 | 12 |
|-|---|---|---|---|---|----|
| NAG MPI   | 340 |  175 |   92 |   62 |   49 |   42 |
| NAG CAF   | 222 |  124 |   79 |   66 |   60 |   64 |
| GNU MPI   | 182 |   95 |   52 |   36 |   31 |   25 |
| GNU CAF   | 183 | 1238 | 3630 | 3824 | 4300 | 4528 |
| Intel MPI | 204 |  107 |   58 |   41 |   33 |   32 |
| Intel CAF | 205 |  DNF |  DNF |  DNF |  DNF |  DNF |

#### Using `NZ=101`

| | 1 | 2 | 4 | 8 |
|-|---|---|---|---|
| Intel MPI | 50 |     35 |   21 |   15 |
| Intel CAF | 49 |  21077 |  27042 |  66040 |
