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

submodule (index_map_type) gather_impl
contains

  module subroutine gather(this, local_data)
    class(index_map), intent(in) :: this
    integer, intent(inout) :: local_data(:)
    call gather_aux(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  subroutine gather_aux(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer, intent(in)  :: onp_data(:)
    integer, intent(out) :: offp_data(:)
    integer :: ierr
    integer, allocatable :: send_buf(:)
    send_buf = onp_data(this%send_index)
    call MPI_Neighbor_alltoallv(send_buf, this%send_counts, this%send_displs, MPI_INTEGER, &
        offp_data, this%recv_counts, this%recv_displs, MPI_INTEGER, this%comm, ierr)
  end subroutine

!  subroutine gather2(this, local_data)
!    class(index_map), intent(in) :: this
!    integer, intent(inout) :: local_data(:,:)
!    call gather2_aux(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
!  end subroutine
!
!  subroutine gather2_aux(this, onp_data, offp_data)
!    class(index_map), intent(in) :: this
!    integer, intent(in)  :: onp_data(:,:)
!    integer, intent(out) :: offp_data(:,:)
!    integer :: ierr
!    integer, allocatable :: send_buf(:,:)
!    integer :: block_type
!    send_buf = onp_data(:,this%send_index)
!    call MPI_Type_contiguous(size(onp_data,dim=1), MPI_INTEGER, block_type, ierr)
!    call MPI_Type_commit(block_type, ierr)
!    call MPI_Neighbor_alltoallv(send_buf, this%send_counts, this%send_displs, block_type, &
!        offp_data, this%recv_counts, this%recv_displs, block_type, this%comm, ierr)
!    call MPI_Type_free(block_type, ierr)
!  end subroutine
!
!  subroutine gather3(this, local_data)
!    class(index_map), intent(in) :: this
!    integer, intent(inout) :: local_data(:,:,:)
!    call gather3_aux(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
!  end subroutine
!
!  subroutine gather3_aux(this, onp_data, offp_data)
!    class(index_map), intent(in) :: this
!    integer, intent(in)  :: onp_data(:,:,:)
!    integer, intent(out) :: offp_data(:,:,:)
!    integer :: ierr
!    integer, allocatable :: send_buf(:,:,:)
!    integer :: block_type
!    send_buf = onp_data(:,:,this%send_index)
!    call MPI_Type_contiguous(size(onp_data(:,:,1)), MPI_INTEGER, block_type, ierr)
!    call MPI_Type_commit(block_type, ierr)
!    call MPI_Neighbor_alltoallv(send_buf, this%send_counts, this%send_displs, block_type, &
!        offp_data, this%recv_counts, this%recv_displs, block_type, this%comm, ierr)
!    call MPI_Type_free(block_type, ierr)
!  end subroutine
!
!  subroutine scatter_sum(this, local_data)
!    class(index_map), intent(in) :: this
!    integer, intent(inout) :: local_data(:)
!    call scatter_sum_aux(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
!  end subroutine
!
!!TODO: %send_*, %recv_* named for gather, but opposite for scatter. Better names?
!!TODO: They are associated with the onp and offp sides, respectively
!
!  subroutine scatter_sum_aux(this, onp_data, offp_data)
!    class(index_map), intent(in) :: this
!    integer, intent(inout) :: onp_data(:)
!    integer, intent(in) :: offp_data(:)
!    integer :: j, ierr
!    integer, allocatable :: onp_buf(:)
!    allocate(onp_buf, mold=this%send_index) !TODO: change name send_index -> onp_index
!    call MPI_Neighbor_alltoallv(offp_data, this%recv_counts, this%recv_displs, MPI_INTEGER, &
!        onp_buf, this%send_counts, this%send_displs, MPI_INTEGER, this%comm, ierr)
!    !NB: This must be done sequentially; SEND_INDEX is a many-to-one mapping.
!    do j = 1, size(onp_buf)
!      onp_data(this%send_index(j)) = onp_data(this%send_index(j)) + onp_buf(j)
!    end do
!  end subroutine
!
!  subroutine scatter_max(this, local_data)
!    class(index_map), intent(in) :: this
!    integer, intent(inout) :: local_data(:)
!    call scatter_max_aux(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
!  end subroutine
!
!!TODO: %send_*, %recv_* named for gather, but opposite for scatter. Better names?
!!TODO: They are associated with the onp and offp sides, respectively
!
!  subroutine scatter_max_aux(this, onp_data, offp_data)
!    class(index_map), intent(in) :: this
!    integer, intent(inout) :: onp_data(:)
!    integer, intent(in) :: offp_data(:)
!    integer :: j, ierr
!    integer, allocatable :: onp_buf(:)
!    allocate(onp_buf, mold=this%send_index) !TODO: change name send_index -> onp_index
!    call MPI_Neighbor_alltoallv(offp_data, this%recv_counts, this%recv_displs, MPI_INTEGER, &
!        onp_buf, this%send_counts, this%send_displs, MPI_INTEGER, this%comm, ierr)
!    !NB: This must be done sequentially; SEND_INDEX is a many-to-one mapping.
!    do j = 1, size(onp_buf)
!      onp_data(this%send_index(j)) = max(onp_data(this%send_index(j)), onp_buf(j))
!    end do
!  end subroutine

end submodule gather_impl
