## Test Configuration
* NAG Fortran 7.1.7106, MPICH 4.0.1
  - MPI version:
    ```
    cmake -DCMAKE_Fortran_FLAGS="-O3 -DNDEBUG" ..
    mpirun -np <N> -bind-to=rr disk-fv-parallel
    ```
  - Coarray version:
    ```
    cmake -DCMAKE_Fortran_FLAGS="-O3 -DNDEBUG" -DUSE_CAF=Yes ..
    NAGFORTRAN_NUM_IMAGES=<N> disk-fv-parallel
    ```
* GFortran 11.2.0, MPICH 4.0.1, OpenCoarrays 2.9.2-13-g235167d
  - MPI version:
    ```
    cmake -DCMAKE_Fortran_FLAGS="-O3 -DNDEBUG" ..
    mpirun -np <N> -bind-to=rr disk-fv-parallel
    ```
  - Coarray version:
    ```
    FC=caf cmake -DCMAKE_Fortran_FLAGS="-O3 -DNDEBUG" -DUSE_CAF=Yes ..
    cafrun -n <N> -bind-to=rr disk-fv-parallel
    ```
* Intel Fortran/MPI 2021.5.0
  - Environment settings for both MPI and CAF:
    ```
    export I_MPI_FABRICS=shm:ofi
    export I_MPI_PIN_PROCESSOR_LIST=0-11
    ```
  - MPI version:
    ```
    cmake -DCMAKE_Fortran_FLAGS="-O3 -DNDEBUG" ..
    mpirun -np <N> disk-fv-parallel
    ```
  - Coarray version:
    ```
    cmake -DCMAKE_Fortran_FLAGS="-O3 -DNDEBUG" -DUSE_CAF=Yes ..
    FOR_COARRAY_NUM_IMAGES=<N> disk-fv-parallel
    ```
