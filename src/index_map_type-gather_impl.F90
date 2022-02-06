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

submodule(index_map_type) gather_impl
contains

  module subroutine gather1_int32_r1(this, local_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(inout) :: local_data(:)
    call gather2_int32_r1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gather2_int32_r1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(in)  :: onp_data(:)
    integer(int32), intent(out) :: offp_data(:)
    integer :: ierr
    integer(int32), allocatable :: send_buf(:)
    send_buf = onp_data(this%send_index)
    call MPI_Neighbor_alltoallv(send_buf, this%send_counts, this%send_displs, MPI_INTEGER4, &
        offp_data, this%recv_counts, this%recv_displs, MPI_INTEGER4, this%comm, ierr)
  end subroutine

  module subroutine gather1_int32_r2(this, local_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(inout) :: local_data(:,:)
    call gather2_int32_r2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gather2_int32_r2(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(in)  :: onp_data(:,:)
    integer(int32), intent(out) :: offp_data(:,:)
    integer :: ierr
    integer(int32), allocatable :: send_buf(:,:)
    integer :: block_type
    send_buf = onp_data(:,this%send_index)
    call MPI_Type_contiguous(size(onp_data,dim=1), MPI_INTEGER4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(send_buf, this%send_counts, this%send_displs, block_type, &
        offp_data, this%recv_counts, this%recv_displs, block_type, this%comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

  module subroutine gather1_int32_r3(this, local_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(inout) :: local_data(:,:,:)
    call gather2_int32_r3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gather2_int32_r3(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    integer(int32), intent(in)  :: onp_data(:,:,:)
    integer(int32), intent(out) :: offp_data(:,:,:)
    integer :: ierr
    integer(int32), allocatable :: send_buf(:,:,:)
    integer :: block_type
    send_buf = onp_data(:,:,this%send_index)
    call MPI_Type_contiguous(size(onp_data(:,:,1)), MPI_INTEGER4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(send_buf, this%send_counts, this%send_displs, block_type, &
        offp_data, this%recv_counts, this%recv_displs, block_type, this%comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

  module subroutine gather1_real64_r1(this, local_data)
    class(index_map), intent(in) :: this
    real(real64), intent(inout) :: local_data(:)
    call gather2_real64_r1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gather2_real64_r1(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(real64), intent(in)  :: onp_data(:)
    real(real64), intent(out) :: offp_data(:)
    integer :: ierr
    real(real64), allocatable :: send_buf(:)
    send_buf = onp_data(this%send_index)
    call MPI_Neighbor_alltoallv(send_buf, this%send_counts, this%send_displs, MPI_REAL8, &
        offp_data, this%recv_counts, this%recv_displs, MPI_REAL8, this%comm, ierr)
  end subroutine

  module subroutine gather1_real64_r2(this, local_data)
    class(index_map), intent(in) :: this
    real(real64), intent(inout) :: local_data(:,:)
    call gather2_real64_r2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gather2_real64_r2(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(real64), intent(in)  :: onp_data(:,:)
    real(real64), intent(out) :: offp_data(:,:)
    integer :: ierr
    real(real64), allocatable :: send_buf(:,:)
    integer :: block_type
    send_buf = onp_data(:,this%send_index)
    call MPI_Type_contiguous(size(onp_data,dim=1), MPI_REAL8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(send_buf, this%send_counts, this%send_displs, block_type, &
        offp_data, this%recv_counts, this%recv_displs, block_type, this%comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

  module subroutine gather1_real64_r3(this, local_data)
    class(index_map), intent(in) :: this
    real(real64), intent(inout) :: local_data(:,:,:)
    call gather2_real64_r3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gather2_real64_r3(this, onp_data, offp_data)
    class(index_map), intent(in) :: this
    real(real64), intent(in)  :: onp_data(:,:,:)
    real(real64), intent(out) :: offp_data(:,:,:)
    integer :: ierr
    real(real64), allocatable :: send_buf(:,:,:)
    integer :: block_type
    send_buf = onp_data(:,:,this%send_index)
    call MPI_Type_contiguous(size(onp_data(:,:,1)), MPI_REAL8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Neighbor_alltoallv(send_buf, this%send_counts, this%send_displs, block_type, &
        offp_data, this%recv_counts, this%recv_displs, block_type, this%comm, ierr)
    call MPI_Type_free(block_type, ierr)
  end subroutine

end submodule gather_impl
