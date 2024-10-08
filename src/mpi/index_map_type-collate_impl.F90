!! Implementation of INDEX_MAP Collate Subroutines
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

submodule(index_map_type) collate_impl
implicit none
contains

!!!! RANK-1 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine coll_i4_1(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: src(:)
    integer(i4), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_INTEGER4, &
        dest, this%counts, this%displs, MPI_INTEGER4, this%root, this%comm, ierr)
  end subroutine

  module subroutine coll_i8_1(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i8), intent(in) :: src(:)
    integer(i8), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_INTEGER8, &
        dest, this%counts, this%displs, MPI_INTEGER8, this%root, this%comm, ierr)
  end subroutine

  module subroutine coll_r4_1(this, src, dest)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: src(:)
    real(r4), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_REAL4, &
        dest, this%counts, this%displs, MPI_REAL4, this%root, this%comm, ierr)
  end subroutine

  module subroutine coll_r8_1(this, src, dest)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: src(:)
    real(r8), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_REAL8, &
        dest, this%counts, this%displs, MPI_REAL8, this%root, this%comm, ierr)
  end subroutine

  module subroutine coll_c4_1(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r4), intent(in) :: src(:)
    complex(r4), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_COMPLEX8, &
        dest, this%counts, this%displs, MPI_COMPLEX8, this%root, this%comm, ierr)
  end subroutine

  module subroutine coll_c8_1(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r8), intent(in) :: src(:)
    complex(r8), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_COMPLEX16, &
        dest, this%counts, this%displs, MPI_COMPLEX16, this%root, this%comm, ierr)
  end subroutine

  module subroutine coll_dl_1(this, src, dest)
    class(index_map), intent(in) :: this
    logical, intent(in) :: src(:)
    logical, intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    call MPI_Gatherv(src, this%onp_size, MPI_LOGICAL, &
        dest, this%counts, this%displs, MPI_LOGICAL, this%root, this%comm, ierr)
  end subroutine

!!!! RANK-2 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! This auxiliary subroutine performs the collate operation for elements that
  !! are blocks of values of the specified MPI data type. The block size need
  !! only be specified on the root process. Note that the array dummy arguments
  !! are assumed-size, so if the actual arguments are not contiguous, copy-in/
  !! copy-out of contiguous temporaries will occur.

  subroutine coll_aux2(this, block_size, mpi_type, src, dest)

    type(index_map), intent(in) :: this
    integer, intent(inout) :: block_size
    integer, intent(in) :: mpi_type
    type(*), intent(in) :: src(*)
    type(*), intent(inout) :: dest(*)

    interface ! explicit interface needed to pass assumed-type arguments
      subroutine MPI_Gatherv(sendbuf, sendcount, sendtype, &
          recvbuf, recvcounts, displs, recvtype, root, comm, ierr)
        type(*), intent(in) :: sendbuf(*)
        type(*), intent(inout) :: recvbuf(*)
        integer, intent(in) :: sendcount, sendtype, recvcounts(*), displs(*), recvtype, root, comm
        integer, intent(out) :: ierr
      end subroutine
    end interface

    integer :: ierr
    integer :: block_type

    call MPI_Bcast(block_size, 1, MPI_INTEGER, this%root, this%comm, ierr)
    call MPI_Type_contiguous(block_size, mpi_type, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(src, this%onp_size, block_type, &
        dest, this%counts, this%displs, block_type, this%root, this%comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine coll_aux2

  module subroutine coll_i4_2(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: src(:,:)
    integer(i4), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call coll_aux2(this, block_size, MPI_INTEGER4, src, dest)
  end subroutine

  module subroutine coll_i8_2(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i8), intent(in) :: src(:,:)
    integer(i8), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call coll_aux2(this, block_size, MPI_INTEGER8, src, dest)
  end subroutine

  module subroutine coll_r4_2(this, src, dest)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: src(:,:)
    real(r4), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call coll_aux2(this, block_size, MPI_REAL4, src, dest)
  end subroutine

  module subroutine coll_r8_2(this, src, dest)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: src(:,:)
    real(r8), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call coll_aux2(this, block_size, MPI_REAL8, src, dest)
  end subroutine

  module subroutine coll_c4_2(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r4), intent(in) :: src(:,:)
    complex(r4), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call coll_aux2(this, block_size, MPI_COMPLEX8, src, dest)
  end subroutine

  module subroutine coll_c8_2(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r8), intent(in) :: src(:,:)
    complex(r8), intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call coll_aux2(this, block_size, MPI_COMPLEX16, src, dest)
  end subroutine

  module subroutine coll_dl_2(this, src, dest)
    class(index_map), intent(in) :: this
    logical, intent(in) :: src(:,:)
    logical, intent(inout) :: dest(:,:)
    integer :: block_size
    ASSERT(size(src,2) >= this%onp_size)
    ASSERT(size(dest,2) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,1))
    call coll_aux2(this, block_size, MPI_LOGICAL, src, dest)
  end subroutine

!!!! RANK-3 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine coll_i4_3(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i4), intent(in) :: src(:,:,:)
    integer(i4), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call coll_aux2(this, block_size, MPI_INTEGER4, src, dest)
  end subroutine

  module subroutine coll_i8_3(this, src, dest)
    class(index_map), intent(in) :: this
    integer(i8), intent(in) :: src(:,:,:)
    integer(i8), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call coll_aux2(this, block_size, MPI_INTEGER8, src, dest)
  end subroutine

  module subroutine coll_r4_3(this, src, dest)
    class(index_map), intent(in) :: this
    real(r4), intent(in) :: src(:,:,:)
    real(r4), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call coll_aux2(this, block_size, MPI_REAL4, src, dest)
  end subroutine

  module subroutine coll_r8_3(this, src, dest)
    class(index_map), intent(in) :: this
    real(r8), intent(in) :: src(:,:,:)
    real(r8), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call coll_aux2(this, block_size, MPI_REAL8, src, dest)
  end subroutine

  module subroutine coll_c4_3(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r4), intent(in) :: src(:,:,:)
    complex(r4), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call coll_aux2(this, block_size, MPI_COMPLEX8, src, dest)
  end subroutine

  module subroutine coll_c8_3(this, src, dest)
    class(index_map), intent(in) :: this
    complex(r8), intent(in) :: src(:,:,:)
    complex(r8), intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call coll_aux2(this, block_size, MPI_COMPLEX16, src, dest)
  end subroutine

  module subroutine coll_dl_3(this, src, dest)
    class(index_map), intent(in) :: this
    logical, intent(in) :: src(:,:,:)
    logical, intent(inout) :: dest(:,:,:)
    integer :: block_size
    ASSERT(size(src,3) >= this%onp_size)
    ASSERT(size(dest,3) >= merge(this%global_size,0,this%is_root))
    if (this%global_size == 0) return ! nothing to do
    if (this%is_root) block_size = size(dest(:,:,1))
    call coll_aux2(this, block_size, MPI_LOGICAL, src, dest)
  end subroutine

end submodule collate_impl
