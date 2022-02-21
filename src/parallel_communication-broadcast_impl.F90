!! Implementation of PARALLEL_COMMUNICATION BROADCAST Procedures
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

submodule(parallel_communication) broadcast_impl
implicit none
contains

!!!! BROADCAST SCALAR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine bcast_i1_0(scalar)
    integer(int8), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine bcast_i4_0(scalar)
    integer(int32), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine bcast_i8_0(scalar)
    integer(int64), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine bcast_r4_0(scalar)
    real(real32), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine bcast_r8_0(scalar)
    real(real64), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine bcast_log_0(scalar)
    logical, intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, 1, MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine bcast_char_0(scalar)
    character(*), intent(inout) :: scalar
    integer :: ierr
    call MPI_Bcast(scalar, len(scalar), MPI_CHARACTER, root, comm, ierr)
  end subroutine

!!!! BROADCAST RANK-1 ARRAY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine bcast_i1_1(vector)
    integer(int8), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine bcast_i4_1(vector)
    integer(int32), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine bcast_i8_1(vector)
    integer(int64), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine bcast_r4_1(vector)
    real(real32), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine bcast_r8_1(vector)
    real(real64), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine bcast_log_1(vector)
    logical, intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine bcast_char_1(vector)
    character(*), intent(inout) :: vector(:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector)*len(vector), MPI_CHARACTER, root, comm, ierr)
  end subroutine

!!!! BROADCAST RANK-2 ARRAY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine bcast_i1_2(vector)
    integer(int8), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine bcast_i4_2(vector)
    integer(int32), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine bcast_i8_2(vector)
    integer(int64), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine bcast_r4_2(vector)
    real(real32), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine bcast_r8_2(vector)
    real(real64), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine bcast_log_2(vector)
    logical, intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine bcast_char_2(vector)
    character(*), intent(inout) :: vector(:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector)*len(vector), MPI_CHARACTER, root, comm, ierr)
  end subroutine

!!!! BROADCAST RANK-3 ARRAY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine bcast_i1_3(vector)
    integer(int8), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER1, root, comm, ierr)
  end subroutine

  module subroutine bcast_i4_3(vector)
    integer(int32), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER4, root, comm, ierr)
  end subroutine

  module subroutine bcast_i8_3(vector)
    integer(int64), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_INTEGER8, root, comm, ierr)
  end subroutine

  module subroutine bcast_r4_3(vector)
    real(real32), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL4, root, comm, ierr)
  end subroutine

  module subroutine bcast_r8_3(vector)
    real(real64), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_REAL8, root, comm, ierr)
  end subroutine

  module subroutine bcast_log_3(vector)
    logical, intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector), MPI_LOGICAL, root, comm, ierr)
  end subroutine

  module subroutine bcast_char_3(vector)
    character(*), intent(inout) :: vector(:,:,:)
    integer :: ierr
    call MPI_Bcast(vector, size(vector)*len(vector), MPI_CHARACTER, root, comm, ierr)
  end subroutine

end submodule broadcast_impl
