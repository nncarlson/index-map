!! Implementation of INDEX_MAP LOCALIZE_INDEX_ARRAY Subroutines
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!! Permission is hereby granted, free of charge, to any person obtaining a
!! copy of this software and associated documentation files (the "Software"),
!! to deal in the Software without restriction, including without limitation
!! the rights to use, copy, modify, merge, publish, distribute, sublicense,
!! and/or sell copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following conditions:
!!
!! The above copyright notice and this permission notice shall be included
!! in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!! DEALINGS IN THE SOFTWARE.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

submodule(index_map_type) localize_impl
use integer_set_type
use integer_map_type
implicit none
contains

  !! Localize a serial rank-1 indexing array G_INDEX given on the root process.
  !! The array is distributed according to the DOMAIN index map and then its
  !! values localized according to the RANGE index map.

  module subroutine localize_index_array_serial_1(domain, g_index, range, l_index, stat)
    class(index_map), intent(inout) :: domain ! inout due to %distribute
    integer, intent(in) :: g_index(:)
    class(index_map), intent(inout) :: range
    integer, allocatable, intent(out) :: l_index(:)
    integer, intent(out), optional :: stat
    if (domain%global_size == 0) then
      allocate(l_index(0))
      if (present(stat)) stat = 0
      return
    end if
    ASSERT(size(g_index) >= merge(domain%global_size,0,domain%is_root))
    allocate(l_index(domain%local_size))
    call domain%distribute(g_index, l_index)
    if (allocated(domain%offp_index)) call domain%gather_offp(l_index)
    call range%localize_index_array(l_index, stat)
  end subroutine

  !! Localize a serial rank-2 indexing array G_INDEX given on the root process.
  !! The array is distributed according to the DOMAIN index map and then its
  !! values localized according to the RANGE index map.

  module subroutine localize_index_array_serial_2(domain, g_index, range, l_index, stat)
    class(index_map), intent(inout) :: domain ! inout due to %distribute
    integer, intent(in) :: g_index(:,:)
    class(index_map), intent(inout) :: range
    integer, allocatable, intent(out) :: l_index(:,:)
    integer, intent(out), optional :: stat
    ASSERT(size(g_index,dim=2) >= merge(domain%global_size,0,domain%is_root))
    if (domain%global_size == 0) then
      allocate(l_index(size(g_index,dim=1),0))  !TODO? Get dim=1 size from root?
      if (present(stat)) stat = 0
      return
    end if
    allocate(l_index(size(g_index,1),domain%local_size))
    call domain%distribute(g_index, l_index)
    if (allocated(domain%offp_index)) call domain%gather_offp(l_index)
    call range%localize_index_array(l_index, stat)
  end subroutine

  !! Localize a distributed rank-1 indexing array INDEX with respect to RANGE.
  !! RANGE will be modified to include off-process references if necessary.
  !! It is an error (stat==-1) if RANGE already includes off-process references
  !! and they do not completely satisfy the off-process references of INDEX.

  module subroutine localize_index_array_dist_1(range, index, stat)

    class(index_map), intent(inout) :: range
    integer, intent(inout) :: index(:)
    integer, intent(out), optional :: stat

    integer :: n, j, ierr
    type(integer_set) :: offp_set
    type(integer_map) :: lid

    ASSERT(minval(index) >= 0)  ! 0 is allowed as a special ID
    ASSERT(maxval(index) <= range%global_size)

    !! Identify all off-process index references with respect to RANGE.
    do j = 1, size(index)
      n = index(j)
      if (n == 0) cycle ! special ID that is preseved
      if (n < range%first_gid) then
        call offp_set%add(n)
      else if (n > range%last_gid) then
        call offp_set%add(n)
      end if
    end do

    !! Adding additional off-process indices to an index_map object that already
    !! has some is not supported. If RANGE has off-process indices, check that its
    !! existing list includes all the ids identified above.
    if (allocated(range%offp_index)) then
      call offp_set%remove(range%offp_index)
      call MPI_Allreduce(offp_set%size(), n, 1, MPI_INTEGER, MPI_MAX, range%comm, ierr)
      if (n > 0) then
        if (present(stat)) then ! caller handles error
          stat = 1
          return
        end if
        block ! we handle error
          use,intrinsic :: iso_fortran_env, only: error_unit
          if (range%is_root) write(error_unit,'(a)') &
              'ERROR: index_map%localize_index_array: extending existing offp_index array'
          call MPI_Finalize(ierr)
          error stop
        end block
      end if
    else
      call MPI_Allreduce(offp_set%size(), n, 1, MPI_INTEGER, MPI_MAX, range%comm, ierr)
      if (n > 0) call add_offp_index_set(range, offp_set)
    end if

    !! Mapping from off-process global ids to local ids with respect to range.
    if (allocated(range%offp_index)) then
      do j = 1, size(range%offp_index)
        call lid%set(key=range%offp_index(j), val=j+range%onp_size)
      end do
    end if

    !! Map index array global ids to their local ids.
    do j = 1, size(index)
      n = index(j)
      if (n == 0) cycle ! special ID that is preseved
      if (n < range%first_gid) then
        index(j) = lid%val(n)
      else if (n > range%last_gid) then
        index(j) = lid%val(n)
      else
        index(j) = n - range%first_gid + 1
      end if
    end do

  end subroutine localize_index_array_dist_1

  !! Localize a distributed rank-2 indexing array INDEX with respect to RANGE.
  !! RANGE will be modified to include off-process references if necessary.
  !! It is an error (stat==-1) if RANGE already includes off-process references
  !! and they do not completely satisfy the off-process references of INDEX.

  module subroutine localize_index_array_dist_2(range, index, stat)
    class(index_map), intent(inout) :: range
    integer, contiguous, intent(inout), target :: index(:,:)
    integer, intent(out), optional :: stat
    integer, pointer :: flat(:)
    flat(1:size(index)) => index
    call range%localize_index_array(flat, stat)
  end subroutine

  module subroutine localize_index_struct_serial(domain, g_count, g_index, range, l_count, l_index, stat)

    class(index_map), intent(inout) :: domain ! inout due to %distribute
    integer, intent(in) :: g_index(:), g_count(:)
    class(index_map), intent(inout) :: range
    integer, allocatable, intent(out) :: l_index(:), l_count(:)
    integer, intent(out), optional :: stat

    integer :: i, j, n, ierr
    type(index_map) :: flat_map
    integer, allocatable :: offset(:), flat_offp(:)

    if (domain%is_root) then
      ASSERT( minval(g_count) >= 0 )
      ASSERT( sum(g_count) == size(g_index) )
      ASSERT( size(g_count) == domain%global_size )
      ASSERT( minval(g_index) >= 0 )
      ASSERT( maxval(g_index) <= range%global_size )
    end if

    if (present(stat)) stat = 0

    if (domain%global_size == 0) then ! empty corner case
      allocate(l_count(0), l_index(0))
      return
    end if

    !! Distribute the global count array according to the domain index map.
    allocate(l_count(domain%local_size))
    call domain%distribute(g_count, l_count)
    if (allocated(domain%offP_index)) call domain%gather_offp(l_count)

    !! To distribute the elements of the ragged g_index array we create a
    !! temporary flat index space map corresponding to elements of g_index.
    call flat_map%init(domain%comm, sum(l_count(1:domain%onp_size)))

    !! If the domain index map includes off-process indices, the corresponding
    !! elements of g_index must also be distributed. We accomplish this with
    !! a gather operation on the flat index space map. But first we need to
    !! generate just what those off-process indices are and add them to the
    !! flat index space map. This is a bit subtle. They originate from the
    !! off-process indices in the domain map, and we do this in parallel
    !! locally on each rank.
    if (allocated(domain%offp_index)) then
      call MPI_Scan(flat_map%onp_size, n, 1, MPI_INTEGER, MPI_SUM, flat_map%comm, ierr)
      n = n - flat_map%onp_size ! exclusive scan of flat_map%onp_size
      allocate(offset(domain%local_size))
      if (size(offset) > 0) then
        offset(1) = n
        do j = 2, domain%onp_size
          offset(j) = offset(j-1) + l_count(j-1)
        end do
        call domain%gather_offp(offset)
      end if
      allocate(flat_offp(sum(l_count(domain%onp_size+1:))))
      n = 0
      do j = domain%onp_size + 1, domain%local_size
        flat_offp(n+1:n+l_count(j)) = [(offset(j)+i, i=1,l_count(j))]
        n = n + l_count(j)
      end do
      call add_offp_index(flat_map, flat_offp)
      deallocate(flat_offp, offset)
    endif

    !! Distribute the global index array according to the flat index space map.
    !! The end result is a ragged g_index distributed according to the domain
    !! index map, including ragged off-process indices, if any.
    allocate(l_index(flat_map%local_size))
    call flat_map%distribute(g_index, l_index)
    if (allocated(flat_map%offp_index)) call flat_map%gather_offp(l_index)

    call range%localize_index_array(l_index, stat)

  end subroutine localize_index_struct_serial

end submodule localize_impl
