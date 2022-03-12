!! PARALLEL_COMMUNICATION
!!
!! This is a reimplementation of the original Truchas module that replaces
!! internal use of PGSLib with direct use of MPI. All interfaces are unchanged.
!!
!! * Retains numbering processes (PEs) starting with 1.
!! * The IO process is PE 1 (MPI rank 0). Was user-choice with PGSLib.
!! * User may initialize MPI; if not, module will initialize. Whoever
!!   initializes MPI is responsible for finalizing MPI.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright 2022  Neil N. Carlson <neil.n.carlson@gmail.com>
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

module parallel_communication

  use mpi
  use,intrinsic :: iso_fortran_env, only: int8, int32, int64, real32, real64
  implicit none
  private

  public :: init_parallel_communication, halt_parallel_communication, abort_parallel_communication
  public :: broadcast, distribute, collate
  public :: global_any, global_all, global_count
  public :: global_sum, global_minval, global_maxval, global_dot_product

  integer, parameter :: root = 0
  integer, parameter, public :: comm = MPI_COMM_WORLD
  integer, parameter, public :: io_pe = root + 1

  integer, public, protected :: npe = 1
  integer, public, protected :: this_pe = io_pe
  logical, public, protected :: is_iop = .true.

  interface broadcast
    module subroutine bcast_i1_0(scalar)
      integer(int8), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_i4_0(scalar)
      integer(int32), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_i8_0(scalar)
      integer(int64), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_r4_0(scalar)
      real(real32), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_r8_0(scalar)
      real(real64), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_log_0(scalar)
      logical, intent(inout) :: scalar
    end subroutine
    module subroutine bcast_char_0(scalar)
      character(*), intent(inout) :: scalar
    end subroutine
    module subroutine bcast_i1_1(vector)
      integer(int8), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_i4_1(vector)
      integer(int32), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_i8_1(vector)
      integer(int64), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_r4_1(vector)
      real(real32), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_r8_1(vector)
      real(real64), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_log_1(vector)
      logical, intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_char_1(vector)
      character(*), intent(inout) :: vector(:)
    end subroutine
    module subroutine bcast_i1_2(vector)
      integer(int8), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_i4_2(vector)
      integer(int32), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_i8_2(vector)
      integer(int64), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_r4_2(vector)
      real(real32), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_r8_2(vector)
      real(real64), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_log_2(vector)
      logical, intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_char_2(vector)
      character(*), intent(inout) :: vector(:,:)
    end subroutine
    module subroutine bcast_i1_3(vector)
      integer(int8), intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_i4_3(vector)
      integer(int32), intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_i8_3(vector)
      integer(int64), intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_r4_3(vector)
      real(real32), intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_r8_3(vector)
      real(real64), intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_log_3(vector)
      logical, intent(inout) :: vector(:,:,:)
    end subroutine
    module subroutine bcast_char_3(vector)
      character(*), intent(inout) :: vector(:,:,:)
    end subroutine
  end interface

  interface distribute
    module subroutine dist_i1_0(src, dest)
      integer(int8), intent(in)  :: src(:)
      integer(int8), intent(out) :: dest
    end subroutine
    module subroutine dist_i4_0(src, dest)
      integer(int32), intent(in)  :: src(:)
      integer(int32), intent(out) :: dest
    end subroutine
    module subroutine dist_i8_0(src, dest)
      integer(int64), intent(in)  :: src(:)
      integer(int64), intent(out) :: dest
    end subroutine
    module subroutine dist_r4_0(src, dest)
      real(real32), intent(in)  :: src(:)
      real(real32), intent(out) :: dest
    end subroutine
    module subroutine dist_r8_0(src, dest)
      real(real64), intent(in)  :: src(:)
      real(real64), intent(out) :: dest
    end subroutine
    module subroutine dist_log_0(src, dest)
      logical, intent(in)  :: src(:)
      logical, intent(out) :: dest
    end subroutine
    module subroutine dist_i1_1(src, dest)
      integer(int8), intent(in)  :: src(:)
      integer(int8), intent(out) :: dest(:)
    end subroutine
    module subroutine dist_i4_1(src, dest)
      integer(int32), intent(in)  :: src(:)
      integer(int32), intent(out) :: dest(:)
    end subroutine
    module subroutine dist_i8_1(src, dest)
      integer(int64), intent(in)  :: src(:)
      integer(int64), intent(out) :: dest(:)
    end subroutine
    module subroutine dist_r4_1(src, dest)
      real(real32), intent(in)  :: src(:)
      real(real32), intent(out) :: dest(:)
    end subroutine
    module subroutine dist_r8_1(src, dest)
      real(real64), intent(in)  :: src(:)
      real(real64), intent(out) :: dest(:)
    end subroutine
    module subroutine dist_log_1(src, dest)
      logical, intent(in)  :: src(:)
      logical, intent(out) :: dest(:)
    end subroutine
    module subroutine dist_i1_2(src, dest)
      integer(int8), intent(in)  :: src(:,:)
      integer(int8), intent(out) :: dest(:,:)
    end subroutine
    module subroutine dist_i4_2(src, dest)
      integer(int32), intent(in)  :: src(:,:)
      integer(int32), intent(out) :: dest(:,:)
    end subroutine
    module subroutine dist_i8_2(src, dest)
      integer(int64), intent(in)  :: src(:,:)
      integer(int64), intent(out) :: dest(:,:)
    end subroutine
    module subroutine dist_r4_2(src, dest)
      real(real32), intent(in)  :: src(:,:)
      real(real32), intent(out) :: dest(:,:)
    end subroutine
    module subroutine dist_r8_2(src, dest)
      real(real64), intent(in)  :: src(:,:)
      real(real64), intent(out) :: dest(:,:)
    end subroutine
    module subroutine dist_log_2(src, dest)
      logical, intent(in)  :: src(:,:)
      logical, intent(out) :: dest(:,:)
    end subroutine
    module subroutine dist_i1_3(src, dest)
      integer(int8), intent(in)  :: src(:,:,:)
      integer(int8), intent(out) :: dest(:,:,:)
    end subroutine
    module subroutine dist_i4_3(src, dest)
      integer(int32), intent(in)  :: src(:,:,:)
      integer(int32), intent(out) :: dest(:,:,:)
    end subroutine
    module subroutine dist_i8_3(src, dest)
      integer(int64), intent(in)  :: src(:,:,:)
      integer(int64), intent(out) :: dest(:,:,:)
    end subroutine
    module subroutine dist_r4_3(src, dest)
      real(real32), intent(in)  :: src(:,:,:)
      real(real32), intent(out) :: dest(:,:,:)
    end subroutine
    module subroutine dist_r8_3(src, dest)
      real(real64), intent(in)  :: src(:,:,:)
      real(real64), intent(out) :: dest(:,:,:)
    end subroutine
    module subroutine dist_log_3(src, dest)
      logical, intent(in)  :: src(:,:,:)
      logical, intent(out) :: dest(:,:,:)
    end subroutine
  end interface

  interface collate
    module subroutine coll_i1_0(src, dest)
      integer(int8), intent(in)  :: src
      integer(int8), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_i4_0(src, dest)
      integer(int32), intent(in)  :: src
      integer(int32), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_i8_0(src, dest)
      integer(int64), intent(in)  :: src
      integer(int64), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_r4_0(src, dest)
      real(real32), intent(in)  :: src
      real(real32), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_r8_0(src, dest)
      real(real64), intent(in)  :: src
      real(real64), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_log_0(src, dest)
      logical, intent(in)  :: src
      logical, intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_char_0(src, dest)
      character(*), intent(in)  :: src
      character(*), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_i1_1(src, dest)
      integer(int8), intent(in)  :: src(:)
      integer(int8), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_i4_1(src, dest)
      integer(int32), intent(in)  :: src(:)
      integer(int32), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_i8_1(src, dest)
      integer(int64), intent(in)  :: src(:)
      integer(int64), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_r4_1(src, dest)
      real(real32), intent(in)  :: src(:)
      real(real32), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_r8_1(src, dest)
      real(real64), intent(in)  :: src(:)
      real(real64), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_log_1(src, dest)
      logical, intent(in)  :: src(:)
      logical, intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_char_1(src, dest)
      character(*), intent(in)  :: src(:)
      character(*), intent(inout) :: dest(:)
    end subroutine
    module subroutine coll_i1_2(src, dest)
      integer(int8), intent(in)  :: src(:,:)
      integer(int8), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine coll_i4_2(src, dest)
      integer(int32), intent(in)  :: src(:,:)
      integer(int32), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine coll_i8_2(src, dest)
      integer(int64), intent(in)  :: src(:,:)
      integer(int64), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine coll_r4_2(src, dest)
      real(real32), intent(in)  :: src(:,:)
      real(real32), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine coll_r8_2(src, dest)
      real(real64), intent(in)  :: src(:,:)
      real(real64), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine coll_log_2(src, dest)
      logical, intent(in)  :: src(:,:)
      logical, intent(inout) :: dest(:,:)
    end subroutine
    module subroutine coll_char_2(src, dest)
      character(*), intent(in)  :: src(:,:)
      character(*), intent(inout) :: dest(:,:)
    end subroutine
  end interface

  interface global_any
    module function any_0(mask) result(global)
      logical, intent(in) :: mask
      logical :: global
    end function
    module function any_1(mask) result(global)
      logical, intent(in) :: mask(:)
      logical :: global
    end function
    module function any_2(mask) result(global)
      logical, intent(in) :: mask(:,:)
      logical :: global
    end function
  end interface

  interface global_all
    module function all_0(mask) result(global)
      logical, intent(in) :: mask
      logical :: global
    end function
    module function all_1(mask) result(global)
      logical, intent(in) :: mask(:)
      logical :: global
    end function
    module function all_2(mask) result(global)
      logical, intent(in) :: mask(:,:)
      logical :: global
    end function
  end interface

  interface global_count
    module function count_0(mask) result(n)
      logical, intent(in) :: mask
      integer :: n
    end function
    module function count_1(mask) result(n)
      logical, intent(in) :: mask(:)
      integer :: n
    end function
    module function count_2(mask) result(n)
      logical, intent(in) :: mask(:,:)
      integer :: n
    end function
  end interface

  interface global_sum
    module function sum_i4_0(a) result(s)
      integer(int32), intent(in) :: a
      integer(int32) :: s
    end function
    module function sum_i8_0(a) result(s)
      integer(int64), intent(in) :: a
      integer(int64) :: s
    end function
    module function sum_r4_0(a) result(s)
      real(real32), intent(in) :: a
      real(real32) :: s
    end function
    module function sum_r8_0(a) result(s)
      real(real64), intent(in) :: a
      real(real64) :: s
    end function
    module function sum_i4_1(a, mask) result(s)
      integer(int32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int32) :: s
    end function
    module function sum_i8_1(a, mask) result(s)
      integer(int64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int64) :: s
    end function
    module function sum_r4_1(a, mask) result(s)
      real(real32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real32) :: s
    end function
    module function sum_r8_1(a, mask) result(s)
      real(real64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real64) :: s
    end function
  end interface

  interface global_minval
    module function minval_i4_0(a) result(v)
      integer(int32), intent(in) :: a
      integer(int32) :: v
    end function
    module function minval_i8_0(a) result(v)
      integer(int64), intent(in) :: a
      integer(int64) :: v
    end function
    module function minval_r4_0(a) result(v)
      real(real32), intent(in) :: a
      real(real32) :: v
    end function
    module function minval_r8_0(a) result(v)
      real(real64), intent(in) :: a
      real(real64) :: v
    end function
    module function minval_i4_1(a, mask) result(v)
      integer(int32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int32) :: v
    end function
    module function minval_i8_1(a, mask) result(v)
      integer(int64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int64) :: v
    end function
    module function minval_r4_1(a, mask) result(v)
      real(real32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real32) :: v
    end function
    module function minval_r8_1(a, mask) result(v)
      real(real64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real64) :: v
    end function
  end interface

  interface global_maxval
    module function maxval_i4_0(a) result(v)
      integer(int32), intent(in) :: a
      integer(int32) :: v
    end function
    module function maxval_i8_0(a) result(v)
      integer(int64), intent(in) :: a
      integer(int64) :: v
    end function
    module function maxval_r4_0(a) result(v)
      real(real32), intent(in) :: a
      real(real32) :: v
    end function
    module function maxval_r8_0(a) result(v)
      real(real64), intent(in) :: a
      real(real64) :: v
    end function
    module function maxval_i4_1(a, mask) result(v)
      integer(int32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int32) :: v
    end function
    module function maxval_i8_1(a, mask) result(v)
      integer(int64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      integer(int64) :: v
    end function
    module function maxval_r4_1(a, mask) result(v)
      real(real32), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real32) :: v
    end function
    module function maxval_r8_1(a, mask) result(v)
      real(real64), intent(in) :: a(:)
      logical, intent(in), optional :: mask(:)
      real(real64) :: v
    end function
  end interface

  interface global_dot_product
    module function dot_prod_r4(a, b) result(dp)
      real(real32), intent(in) :: a(:), b(:)
      real(real32) :: dp
    end function
    module function dot_prod_r8(a, b) result(dp)
      real(real64), intent(in) :: a(:), b(:)
      real(real64) :: dp
    end function
  end interface

  logical :: initialized = .false., flag = .false.

contains

  subroutine init_parallel_communication
    integer :: ierr
    if (initialized) return ! should only be called once
    initialized = .true.
    call MPI_Initialized(flag, ierr)
    if (.not.flag) call MPI_Init(ierr)
    call MPI_Comm_size(comm, npe, ierr)
    call MPI_Comm_rank(comm, this_pe, ierr)
    this_pe = this_pe + 1 ! start numbering at 1
    is_iop = (this_pe == io_pe)
  end subroutine

  subroutine halt_parallel_communication
    integer :: ierr
    if (.not.flag) call MPI_Finalize(ierr)
  end subroutine

  subroutine abort_parallel_communication
    integer :: ierr
    call MPI_Abort(comm, 1, ierr)
  end subroutine

end module parallel_communication
