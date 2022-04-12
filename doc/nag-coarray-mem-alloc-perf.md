### An observation about NAG memory allocation in coarray mode

* Using the **serial** code `disk-fv-serial.F90` (NZ=257)
* Tests using nagfor build 7106 (same results with 7103)
* The `COMPACT_UPDATE` option uses this 1-line update statement
  ```fortran
  u(j) = u_prev(j) + c*(sum(u_prev(cnhbr(:,j))) - 4*u_prev(j))
  ```
  instead of the (now) default explicit do loop:
  ```fortran
  tmp = -4*u_prev(j)
  do k = 1, size(cnhbr,1)
    tmp = tmp + u_prev(cnhbr(k,j))
  end do
  u(j) = u_prev(j) + c*tmp
  ```
* Run with `NAGFORTRAN_NUM_IMAGES=1`

Here are sample runtimes (µsec per step) when compiled with/without `-coarray`
and with/without `-DCOMPACT_UPDATE`.

| -O2 optimization         | without -coarray | with -coarray |
|--------------------------|------------------|---------------|
| without -DCOMPACT_UPDATE | 143              | 86            |
|    with -DCOMPACT_UPDATE | 645              | 1685          |

| -O3 optimization         | without -coarray | with -coarray |
|--------------------------|------------------|---------------|
| without -DCOMPACT_UPDATE | 131              | 73            |
|    with -DCOMPACT_UPDATE | 839              | 1683          |

* Without `-coarray`, it had already been observed that the compact code is
  far slower (4x) than the equivalent do loop version; unlike gfortran and
  Intel, the NAG compiler is not able to optimize the compact code.

* The code makes no use whatsoever of coarrays, so why adding the `-coarray`
  option makes any difference, let alone such a huge difference, is surprising.
  The timings are for just the computational loop; program startup and shutdown
  are not included. The likely explanation is that with `-coarray` a different
  memory allocator is being used.

* But what is bizarre, is that adding `-coarray` greatly reduces the time for
  the do loop version, and greatly increases the time for the compact code.

  In the latter case where the compiler is unable to optimize the array syntax,
  there are temporary allocations/deallocations occurring within the timed
  loop. The great increase of time here points to a much more costly memory
  allocator in `-coarray` mode.

  The former case where there are no allocations/deallocations in the timed
  loop is harder to explain. One possible explanation is that the `-coarray`
  memory allocator produced better arrays with respect to cache conflicts,
  page alignment, etc. I have no data to suggest this would be more than just
  happenstance in this specific case.
