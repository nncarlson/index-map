!! Implementation of core INDEX_MAP_TYPE helper procedures
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

#include "f90_assert.fpp"

submodule(index_map_type) core_impl
implicit none
contains

  module subroutine add_offp_index_set(this, offp_set)

    use integer_set_type

    class(index_map), intent(inout), target :: this
    type(integer_set), intent(inout) :: offp_set

    integer :: i, j, j1, k, n
    integer, allocatable :: offp_image(:), onp_counts(:)
    integer, allocatable :: last[:], offp_counts(:)[:], offsets(:)[:]

    type box
      integer, pointer :: data(:)
    end type
    type(box), allocatable :: buffer[:]

    !TODO? Allow extending an existing %offp_index
    INSIST(.not.allocated(this%offp_index))

    this%offp_index = offp_set
    this%offp_size  = size(this%offp_index)
    this%local_size = this%onp_size + this%offp_size

    !! Determine the image that owns each off-process index (OFFP_IMAGE).
    !! NB: This assumes OFFP_INDEX is ordered; OFFP_IMAGE will be ordered.
    allocate(last[*], offp_image(this%offp_size))
    last = this%last_gid
    sync all
    i = 1
    do j = 1, size(this%offp_index)
      do while (this%offp_index(j) > last[i])
        i = i + 1
        INSIST(i <= this%nproc)
      end do
      INSIST(i /= this_image())
      offp_image(j) = i
    end do

    !! Get the number of off-process indices owned by each image (OFFP_COUNTS).
    !! NB: This assumes the OFFP_IMAGE array is ordered.
    allocate(offp_counts(this%nproc)[*])
    j1 = 1
    do i = 1, this%nproc
      do j = j1, size(offp_image)
        if (offp_image(j) > i) exit
      end do
      offp_counts(i) = j - j1
      j1 = j
    end do

    !! Communicate the count of off-process indices owned by an image to that
    !! image. The result ONP_COUNTS is the number of on-process indices that
    !! are included as off-process indices on each of the other images.
    sync all
    allocate(onp_counts(this%nproc))
    do i = 1, this%nproc
      onp_counts(i) = offp_counts(this_image())[i]
    end do

    !! Offsets into the ONP_INDEX array for each of the image blocks.
    allocate(offsets(this%nproc)[*])
    offsets(1) = 0
    do i = 2, this%nproc
      offsets(i) = offsets(i-1) + onp_counts(i-1)
    end do

    !! Compress the OFFP_COUNTS array, dropping elements with 0 count.
    sync all
    n = count(offp_counts > 0)
    allocate(this%offp_count(n), this%offp_image(n), this%onp_offset(n))
    n = 0
    do i = 1, this%nproc
      if (offp_counts(i) == 0) cycle
      n = n + 1
      this%offp_image(n) = i
      this%offp_count(n) = offp_counts(i)
      this%onp_offset(n) = offsets(this_image())[i]
    end do

    !! Offsets into the OFFP_INDEX array for each of the image blocks.
    sync all
    offsets(1) = 0
    do i = 2, this%nproc
      offsets(i) = offsets(i-1) + offp_counts(i-1)
    end do

    !! Compress the ONP_COUNTS array, dropping elements with 0 count.
    sync all
    n = count(onp_counts > 0)
    allocate(this%onp_count(n), this%onp_image(n), this%offp_offset(n))
    n = 0
    do i = 1, this%nproc
      if (onp_counts(i) == 0) cycle
      n = n + 1
      this%onp_image(n) = i
      this%onp_count(n) = onp_counts(i)
      this%offp_offset(n) = offsets(this_image())[i]
    end do

    !! Communicate the global off-process indices to their owning images
    !! and map to the corresponding local on-process indices (ONP_INDEX).
    n = sum(this%onp_count)
    allocate(this%onp_index(n), buffer[*])
    buffer%data => this%offp_index
    sync all
    n = 0
    do j = 1, size(this%onp_count)
      associate (i => this%onp_image(j), offset => this%offp_offset(j))
        do k = 1, this%onp_count(j)
          n = n + 1
          this%onp_index(n) = buffer[i]%data(offset+k) - this%first_gid + 1
        end do
      end associate
    end do
    ASSERT(all(this%onp_index >= 1 .and. this%onp_index <= this%onp_size))

    ASSERT(gather_offp_verified(this))

  end subroutine add_offp_index_set

end submodule core_impl
