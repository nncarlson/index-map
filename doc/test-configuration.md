## Test Configuration
* NAG Fortran 7.2.7200, MPICH 4.1.2
  - MPI version:
    ```
    cmake .. -DCMAKE_Fortran_FLAGS="-O3"
    mpirun -np <N> -bind-to=rr disk-fv-parallel
    ```
  - Coarray version:
    ```
    cmake .. -DCMAKE_Fortran_FLAGS="-O3" -DUSE_CAF=Yes
    NAGFORTRAN_NUM_IMAGES=<N> disk-fv-parallel
    ```
* GFortran 12.3.0, MPICH 4.0.3, OpenCoarrays 2.10.2
  - MPI version:
    ```
    cmake .. -DCMAKE_Fortran_FLAGS="-O3"
    mpirun -np <N> -bind-to=rr disk-fv-parallel
    ```
  - Coarray version:
    ```
    FC=caf cmake .. -DCMAKE_Fortran_FLAGS="-O3" -DUSE_CAF=Yes
    cafrun -n <N> -bind-to=rr disk-fv-parallel
    ```
* Intel OneAPI 2024.1 (ifort or ifx with Intel MPI)
  - Environment settings for both MPI and CAF:
    ```
    export I_MPI_FABRICS=shm:ofi
    export I_MPI_PIN_PROCESSOR_LIST=0-11
    ```
  - MPI version:
    ```
    cmake .. -DCMAKE_Fortran_FLAGS="-O3"
    mpirun -np <N> disk-fv-parallel
    ```
  - Coarray version:
    ```
    cmake .. -DCMAKE_Fortran_FLAGS="-O3" -DUSE_CAF=Yes
    FOR_COARRAY_NUM_IMAGES=<N> disk-fv-parallel
    ```
* Intel OneAPI 2024.1 (ifort or ifx), MPICH 4.1.2
  - MPI version:
    ```
    cmake .. -DCMAKE_Fortran_FLAGS="-O3"
    mpirun -np <N> -bind-to=rr disk-fv-parallel
