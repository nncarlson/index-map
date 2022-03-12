!! Implementation of PARALLEL_COMMUNICATION COLLATE Procedures
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

submodule(parallel_communication) collate_impl
implicit none
contains

!!!! COLLATE SCALAR DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine coll_i1_0(src, dest)
    integer(int8), intent(in) :: src
    integer(int8), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_INTEGER1, dest, 1, MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine coll_i4_0(src, dest)
    integer(int32), intent(in) :: src
    integer(int32), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_INTEGER4, dest, 1, MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine coll_i8_0(src, dest)
    integer(int64), intent(in) :: src
    integer(int64), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_INTEGER8, dest, 1, MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine coll_r4_0(src, dest)
    real(real32), intent(in) :: src
    real(real32), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_REAL4, dest, 1, MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine coll_r8_0(src, dest)
    real(real64), intent(in) :: src
    real(real64), intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_REAL8, dest, 1, MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine coll_log_0(src, dest)
    logical, intent(in) :: src
    logical, intent(inout) :: dest(:)
    integer :: ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    call MPI_Gather(src, 1, MPI_LOGICAL, dest, 1, MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine coll_char_0(src, dest)
    character(*), intent(in) :: src
    character(*), intent(inout) :: dest(:)
    integer :: string_len, ierr
    ASSERT(size(dest) >= merge(npe,0,is_iop))
    string_len = len(src) ! must be same for all character args on all processes!
    call MPI_Gather(src, string_len, MPI_CHARACTER, dest, string_len, MPI_CHARACTER, root, comm, ierr)
  end subroutine

!!!! COLLATE RANK-1 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !! This auxiliary subroutine generates the COUNTS and DISPLS arrays for the
  !! receive buffer on the IO process, given the size of the send buffer INLEN
  !! on each of the processes. OUTLEN is the required size of the receive buffer.

  subroutine make_counts_displs(inlen, counts, displs, outlen)

    integer, intent(in) :: inlen
    integer, allocatable :: counts(:), displs(:)
    integer, intent(out) :: outlen

    integer :: ierr, j

    allocate(counts(merge(npe,0,is_iop)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)
    outlen = sum(counts)
    call MPI_Bcast(outlen, 1, MPI_INTEGER, root, comm, ierr)
    if (outlen == 0) return ! nothing to do

    allocate(displs(merge(npe,0,is_iop)))
    if (is_iop) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

  end subroutine

  module subroutine coll_i1_1(src, dest)

    integer(int8), intent(in) :: src(:)
    integer(int8), intent(inout) :: dest(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(src, dim=1)
    call make_counts_displs(inlen, counts, displs, outlen)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(dest,dim=1) >= merge(outlen,0,is_iop))

    call MPI_Gatherv(src, inlen, MPI_INTEGER1, &
        dest, counts, displs, MPI_INTEGER1, root, comm, ierr)

  end subroutine

  module subroutine coll_i4_1(src, dest)

    integer(int32), intent(in) :: src(:)
    integer(int32), intent(inout) :: dest(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(src, dim=1)
    call make_counts_displs(inlen, counts, displs, outlen)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(dest,dim=1) >= merge(outlen,0,is_iop))

    call MPI_Gatherv(src, inlen, MPI_INTEGER4, &
        dest, counts, displs, MPI_INTEGER4, root, comm, ierr)

  end subroutine

  module subroutine coll_i8_1(src, dest)

    integer(int64), intent(in) :: src(:)
    integer(int64), intent(inout) :: dest(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(src, dim=1)
    call make_counts_displs(inlen, counts, displs, outlen)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(dest,dim=1) >= merge(outlen,0,is_iop))

    call MPI_Gatherv(src, inlen, MPI_INTEGER8, &
        dest, counts, displs, MPI_INTEGER8, root, comm, ierr)

  end subroutine

  module subroutine coll_r4_1(src, dest)

    real(real32), intent(in) :: src(:)
    real(real32), intent(inout) :: dest(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(src, dim=1)
    call make_counts_displs(inlen, counts, displs, outlen)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(dest,dim=1) >= merge(outlen,0,is_iop))

    call MPI_Gatherv(src, inlen, MPI_REAL4, &
        dest, counts, displs, MPI_REAL4, root, comm, ierr)

  end subroutine

  module subroutine coll_r8_1(src, dest)

    real(real64), intent(in) :: src(:)
    real(real64), intent(inout) :: dest(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(src, dim=1)
    call make_counts_displs(inlen, counts, displs, outlen)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(dest,dim=1) >= merge(outlen,0,is_iop))

    call MPI_Gatherv(src, inlen, MPI_REAL8, &
        dest, counts, displs, MPI_REAL8, root, comm, ierr)

  end subroutine

  module subroutine coll_log_1(src, dest)

    logical, intent(in) :: src(:)
    logical, intent(inout) :: dest(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(src, dim=1)
    call make_counts_displs(inlen, counts, displs, outlen)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(dest,dim=1) >= merge(outlen,0,is_iop))

    call MPI_Gatherv(src, inlen, MPI_LOGICAL, &
        dest, counts, displs, MPI_LOGICAL, root, comm, ierr)

  end subroutine

  module subroutine coll_char_1(src, dest)

    character(*), intent(in) :: src(:)
    character(*), intent(inout) :: dest(:)

    integer :: inlen, outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(src, dim=1)
    call make_counts_displs(inlen, counts, displs, outlen)
    if (outlen == 0) return ! nothing to do

    ASSERT(size(dest,dim=1) >= merge(outlen,0,is_iop))

    call MPI_Type_contiguous(len(src), MPI_CHARACTER, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(src, inlen, block_type, &
        dest, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

!!!! COLLATE RANK-2 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine coll_aux(inlen, destlen, block_size, mpi_type, src, dest)

    integer, intent(in) :: inlen, destlen, block_size
    integer, intent(in) :: mpi_type
    type(*), intent(in) :: src(*)
    type(*), intent(inout) :: dest(*)

    integer :: outlen, ierr, j
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    call make_counts_displs(inlen, counts, displs, outlen)
    if (outlen == 0) return ! nothing to do
    ASSERT(destlen >= merge(outlen,0,is_iop))

    call MPI_Bcast(block_size, 1, MPI_INTEGER, root, comm, ierr)
    call MPI_Type_contiguous(block_size, mpi_type, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(src, inlen, block_type, dest, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_i1_2(src, dest)

    integer(int8), intent(in) :: src(:,:)
    integer(int8), intent(inout) :: dest(:,:)

    integer :: srclen, destlen, block_size

    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call coll_aux(srclen, destlen, block_size, MPI_INTEGER1, src, dest)

  end subroutine

  module subroutine coll_i4_2(src, dest)

    integer(int32), intent(in) :: src(:,:)
    integer(int32), intent(inout) :: dest(:,:)

    integer :: srclen, destlen, block_size

    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call coll_aux(srclen, destlen, block_size, MPI_INTEGER4, src, dest)

  end subroutine

  module subroutine coll_i8_2(src, dest)

    integer(int64), intent(in) :: src(:,:)
    integer(int64), intent(inout) :: dest(:,:)

    integer :: srclen, destlen, block_size

    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call coll_aux(srclen, destlen, block_size, MPI_INTEGER8, src, dest)

  end subroutine

  module subroutine coll_r4_2(src, dest)

    real(real32), intent(in) :: src(:,:)
    real(real32), intent(inout) :: dest(:,:)

    integer :: srclen, destlen, block_size

    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call coll_aux(srclen, destlen, block_size, MPI_REAL4, src, dest)

  end subroutine

  module subroutine coll_r8_2(src, dest)

    real(real64), intent(in) :: src(:,:)
    real(real64), intent(inout) :: dest(:,:)

    integer :: srclen, destlen, block_size

    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call coll_aux(srclen, destlen, block_size, MPI_REAL8, src, dest)

  end subroutine

  module subroutine coll_log_2(src, dest)

    logical, intent(in) :: src(:,:)
    logical, intent(inout) :: dest(:,:)

    integer :: srclen, destlen, block_size

    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = size(dest(:,1))
    call coll_aux(srclen, destlen, block_size, MPI_LOGICAL, src, dest)

  end subroutine

  module subroutine coll_char_2(src, dest)

    character(*), intent(in) :: src(:,:)
    character(*), intent(inout) :: dest(:,:)

    integer :: srclen, destlen, block_size

    srclen = size(src, dim=2)
    destlen = size(dest, dim=2)
    if (is_iop) block_size = len(dest)*size(dest(:,1))
    call coll_aux(srclen, destlen, block_size, MPI_CHARACTER, src, dest)

  end subroutine

end submodule collate_impl
