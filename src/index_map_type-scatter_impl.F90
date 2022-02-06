!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2022  Neil N. Carlson
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

submodule(index_map_type) scatter_impl
contains

  module subroutine scatter1_sum_int32_r1(this, local_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(inout) :: local_data(:)
    call scatter2_sum_int32_r1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

!TODO: %send_*, %recv_* named for gather, but opposite for scatter. Better names?
!TODO: They are associated with the onp and offp sides, respectively

  module subroutine scatter2_sum_int32_r1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(inout) :: onp_data(:)
    integer(int32), intent(in) :: offp_data(:)
    integer :: j, ierr
    integer(int32), allocatable :: onp_buf(:)
    allocate(onp_buf(size(this%send_index))) !TODO: change name send_index -> onp_index
    call MPI_Neighbor_alltoallv(offp_data, this%recv_counts, this%recv_displs, MPI_INTEGER4, &
        onp_buf, this%send_counts, this%send_displs, MPI_INTEGER4, this%comm, ierr)
    !NB: This must be done sequentially; SEND_INDEX is a many-to-one mapping.
    do j = 1, size(onp_buf)
      onp_data(this%send_index(j)) = onp_data(this%send_index(j)) + onp_buf(j)
    end do
  end subroutine

  module subroutine scatter1_min_int32_r1(this, local_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(inout) :: local_data(:)
    call scatter2_min_int32_r1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

!TODO: %send_*, %recv_* named for gather, but opposite for scatter. Better names?
!TODO: They are associated with the onp and offp sides, respectively

  module subroutine scatter2_min_int32_r1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(inout) :: onp_data(:)
    integer(int32), intent(in) :: offp_data(:)
    integer :: j, ierr
    integer(int32), allocatable :: onp_buf(:)
    allocate(onp_buf(size(this%send_index))) !TODO: change name send_index -> onp_index
    call MPI_Neighbor_alltoallv(offp_data, this%recv_counts, this%recv_displs, MPI_INTEGER4, &
        onp_buf, this%send_counts, this%send_displs, MPI_INTEGER4, this%comm, ierr)
    !NB: This must be done sequentially; SEND_INDEX is a many-to-one mapping.
    do j = 1, size(onp_buf)
      onp_data(this%send_index(j)) = min(onp_data(this%send_index(j)), onp_buf(j))
    end do
  end subroutine

  module subroutine scatter1_max_int32_r1(this, local_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(inout) :: local_data(:)
    call scatter2_max_int32_r1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

!TODO: %send_*, %recv_* named for gather, but opposite for scatter. Better names?
!TODO: They are associated with the onp and offp sides, respectively

  module subroutine scatter2_max_int32_r1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(inout) :: onp_data(:)
    integer(int32), intent(in) :: offp_data(:)
    integer :: j, ierr
    integer(int32), allocatable :: onp_buf(:)
    allocate(onp_buf(size(this%send_index))) !TODO: change name send_index -> onp_index
    call MPI_Neighbor_alltoallv(offp_data, this%recv_counts, this%recv_displs, MPI_INTEGER4, &
        onp_buf, this%send_counts, this%send_displs, MPI_INTEGER4, this%comm, ierr)
    !NB: This must be done sequentially; SEND_INDEX is a many-to-one mapping.
    do j = 1, size(onp_buf)
      onp_data(this%send_index(j)) = max(onp_data(this%send_index(j)), onp_buf(j))
    end do
  end subroutine

  module subroutine scatter1_or_l_r1(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:)
    call scatter2_or_l_r1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

!TODO: %send_*, %recv_* named for gather, but opposite for scatter. Better names?
!TODO: They are associated with the onp and offp sides, respectively

  module subroutine scatter2_or_l_r1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: onp_data(:)
    logical, intent(in) :: offp_data(:)
    integer :: j, ierr
    logical, allocatable :: onp_buf(:)
    allocate(onp_buf(size(this%send_index))) !TODO: change name send_index -> onp_index
    call MPI_Neighbor_alltoallv(offp_data, this%recv_counts, this%recv_displs, MPI_LOGICAL, &
        onp_buf, this%send_counts, this%send_displs, MPI_LOGICAL, this%comm, ierr)
    !NB: This must be done sequentially; SEND_INDEX is a many-to-one mapping.
    do j = 1, size(onp_buf)
      onp_data(this%send_index(j)) = onp_data(this%send_index(j)) .or. onp_buf(j)
    end do
  end subroutine

  module subroutine scatter1_and_l_r1(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:)
    call scatter2_and_l_r1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

!TODO: %send_*, %recv_* named for gather, but opposite for scatter. Better names?
!TODO: They are associated with the onp and offp sides, respectively

  module subroutine scatter2_and_l_r1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: onp_data(:)
    logical, intent(in) :: offp_data(:)
    integer :: j, ierr
    logical, allocatable :: onp_buf(:)
    allocate(onp_buf(size(this%send_index))) !TODO: change name send_index -> onp_index
    call MPI_Neighbor_alltoallv(offp_data, this%recv_counts, this%recv_displs, MPI_LOGICAL, &
        onp_buf, this%send_counts, this%send_displs, MPI_LOGICAL, this%comm, ierr)
    !NB: This must be done sequentially; SEND_INDEX is a many-to-one mapping.
    do j = 1, size(onp_buf)
      onp_data(this%send_index(j)) = onp_data(this%send_index(j)) .and. onp_buf(j)
    end do
  end subroutine

end submodule scatter_impl
