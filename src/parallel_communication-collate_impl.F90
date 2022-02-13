!! Implementation of PARALLEL_COMMUNICATION COLLATE Procedures
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!! Use subject to the MIT license: https://opensource.org/licenses/MIT
!!

#include "f90_assert.fpp"

submodule(parallel_communication) collate_impl
implicit none
contains

!!!! COLLATE SCALAR DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine coll_i1_0(scalarv_out, scalar_in)
    integer(int8), intent(out) :: scalarv_out(:)
    integer(int8), intent(in)  :: scalar_in
    integer :: ierr
#ifndef NDEBUG
    if (is_iop) then
      ASSERT(size(scalarv_out) == npe)
    end if
#endif
    call MPI_Gather(scalar_in, 1, MPI_INTEGER1, scalarv_out, 1, MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine coll_i4_0(scalarv_out, scalar_in)
    integer(int32), intent(out) :: scalarv_out(:)
    integer(int32), intent(in)  :: scalar_in
    integer :: ierr
#ifndef NDEBUG
    if (is_iop) then
      ASSERT(size(scalarv_out) == npe)
    end if
#endif
    call MPI_Gather(scalar_in, 1, MPI_INTEGER4, scalarv_out, 1, MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine coll_i8_0(scalarv_out, scalar_in)
    integer(int64), intent(out) :: scalarv_out(:)
    integer(int64), intent(in)  :: scalar_in
    integer :: ierr
#ifndef NDEBUG
    if (is_iop) then
      ASSERT(size(scalarv_out) == npe)
    end if
#endif
    call MPI_Gather(scalar_in, 1, MPI_INTEGER8, scalarv_out, 1, MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine coll_r4_0(scalarv_out, scalar_in)
    real(real32), intent(out) :: scalarv_out(:)
    real(real32), intent(in)  :: scalar_in
    integer :: ierr
#ifndef NDEBUG
    if (is_iop) then
      ASSERT(size(scalarv_out) == npe)
    end if
#endif
    call MPI_Gather(scalar_in, 1, MPI_REAL4, scalarv_out, 1, MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine coll_r8_0(scalarv_out, scalar_in)
    real(real64), intent(out) :: scalarv_out(:)
    real(real64), intent(in)  :: scalar_in
    integer :: ierr
#ifndef NDEBUG
    if (is_iop) then
      ASSERT(size(scalarv_out) == npe)
    end if
#endif
    call MPI_Gather(scalar_in, 1, MPI_REAL8, scalarv_out, 1, MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine coll_log_0(scalarv_out, scalar_in)
    logical, intent(out) :: scalarv_out(:)
    logical, intent(in)  :: scalar_in
    integer :: ierr
#ifndef NDEBUG
    if (is_iop) then
      ASSERT(size(scalarv_out) == npe)
    end if
#endif
    call MPI_Gather(scalar_in, 1, MPI_LOGICAL, scalarv_out, 1, MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine coll_char_0(scalarv_out, scalar_in)
    character(*), intent(out) :: scalarv_out(:)
    character(*), intent(in)  :: scalar_in
    integer :: string_len, ierr
#ifndef NDEBUG
    if (is_iop) then
      ASSERT(size(scalarv_out) == npe)
    end if
#endif
    string_len = len(scalar_in) ! must be same for all character args on all processes!
    call MPI_Gather(scalar_in, string_len, MPI_CHARACTER, scalarv_out, string_len, MPI_CHARACTER, root, comm, ierr)
  end subroutine

!!!! COLLATE RANK-1 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!NB: For some reason PGSLib did not implement an optional lengths argument as
!    for the gather subroutines. Thus here it always is computed internally.

  module subroutine coll_i1_1(vector_out, vector_in)

    integer(int8), intent(out) :: vector_out(:)
    integer(int8), intent(in)  :: vector_in(:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=1) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_INTEGER1, &
        vector_out, counts, displs, MPI_INTEGER1, root, comm, ierr)

  end subroutine

  module subroutine coll_i4_1(vector_out, vector_in)

    integer(int32), intent(out) :: vector_out(:)
    integer(int32), intent(in)  :: vector_in(:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=1) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_INTEGER4, &
        vector_out, counts, displs, MPI_INTEGER4, root, comm, ierr)

  end subroutine

  module subroutine coll_i8_1(vector_out, vector_in)

    integer(int64), intent(out) :: vector_out(:)
    integer(int64), intent(in)  :: vector_in(:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=1) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_INTEGER8, &
        vector_out, counts, displs, MPI_INTEGER8, root, comm, ierr)

  end subroutine

  module subroutine coll_r4_1(vector_out, vector_in)

    real(real32), intent(out) :: vector_out(:)
    real(real32), intent(in)  :: vector_in(:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=1) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_REAL4, &
        vector_out, counts, displs, MPI_REAL4, root, comm, ierr)

  end subroutine

  module subroutine coll_r8_1(vector_out, vector_in)

    real(real64), intent(out) :: vector_out(:)
    real(real64), intent(in)  :: vector_in(:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=1) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_REAL8, &
        vector_out, counts, displs, MPI_REAL8, root, comm, ierr)

  end subroutine

  module subroutine coll_log_1(vector_out, vector_in)

    logical, intent(out) :: vector_out(:)
    logical, intent(in)  :: vector_in(:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=1) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Gatherv(vector_in, inlen, MPI_LOGICAL, &
        vector_out, counts, displs, MPI_LOGICAL, root, comm, ierr)

  end subroutine

  module subroutine coll_char_1(vector_out, vector_in)

    character(*), intent(out) :: vector_out(:)
    character(*), intent(in)  :: vector_in(:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=1)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=1) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Type_contiguous(len(vector_in), MPI_CHARACTER, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

!!!! COLLATE RANK-2 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine coll_i1_2(vector_out, vector_in)

    integer(int8), intent(out) :: vector_out(:,:)
    integer(int8), intent(in)  :: vector_in(:,:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=2) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Type_contiguous(size(vector_in(:,1)), MPI_INTEGER1, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_i4_2(vector_out, vector_in)

    integer(int32), intent(out) :: vector_out(:,:)
    integer(int32), intent(in)  :: vector_in(:,:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=2) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Type_contiguous(size(vector_in(:,1)), MPI_INTEGER4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_i8_2(vector_out, vector_in)

    integer(int64), intent(out) :: vector_out(:,:)
    integer(int64), intent(in)  :: vector_in(:,:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=2) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Type_contiguous(size(vector_in(:,1)), MPI_INTEGER8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_r4_2(vector_out, vector_in)

    real(real32), intent(out) :: vector_out(:,:)
    real(real32), intent(in)  :: vector_in(:,:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=2) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Type_contiguous(size(vector_in(:,1)), MPI_REAL4, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_r8_2(vector_out, vector_in)

    real(real64), intent(out) :: vector_out(:,:)
    real(real64), intent(in)  :: vector_in(:,:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=2) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Type_contiguous(size(vector_in(:,1)), MPI_REAL8, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_log_2(vector_out, vector_in)

    logical, intent(out) :: vector_out(:,:)
    logical, intent(in)  :: vector_in(:,:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=2) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Type_contiguous(size(vector_in(:,1)), MPI_LOGICAL, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

  module subroutine coll_char_2(vector_out, vector_in)

    character(*), intent(out) :: vector_out(:,:)
    character(*), intent(in)  :: vector_in(:,:)

    integer :: inlen, ierr, j
    integer, allocatable :: counts(:), displs(:)
    integer :: block_type

    inlen = size(vector_in, dim=2)

    allocate(counts(merge(npe,0,is_IOP)))
    call MPI_Gather(inlen, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, root, comm, ierr)

#ifndef NDEBUG
    if (is_IOP) then
      ASSERT(size(vector_out,dim=2) >= sum(counts))
    end if
#endif

    allocate(displs(merge(npe,0,is_IOP)))
    if (is_IOP) then
      displs(1) = 0
      do j = 2, npe
        displs(j) = displs(j-1) + counts(j-1)
      end do
    end if

    call MPI_Type_contiguous(len(vector_in)*size(vector_in,1), MPI_CHARACTER, block_type, ierr)
    call MPI_Type_commit(block_type, ierr)
    call MPI_Gatherv(vector_in, inlen, block_type, &
        vector_out, counts, displs, block_type, root, comm, ierr)
    call MPI_Type_free(block_type, ierr)

  end subroutine

end submodule collate_impl
