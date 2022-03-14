!! Unit Tests for PARALLEL_COMMUNICATION Distribute Procedures
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

program main

  use mpi
  use parallel_communication
  use,intrinsic :: iso_fortran_env
  implicit none

  integer :: ierr, status, global_status
  logical :: pass ! local to tests

  call MPI_Init(ierr)
  call init_parallel_communication

  status = 0
  call dist_scalar
  call dist_rank1
  call dist_rank1_zero
  call dist_rank2
  call dist_rank2_zero
  call dist_array_section
  call dist_rank3
  call dist_rank3_zero
  call dist_log_scalar
  call dist_log_rank1
  call dist_log_rank2
  call dist_log_rank3

  call MPI_Allreduce(status, global_status, 1, MPI_INTEGER, MPI_MAX, comm, ierr)

  call MPI_Finalize(ierr)

  if (global_status /= 0) stop 1

contains

  subroutine write_result(pass, name)
    logical, intent(in) :: pass
    character(*), intent(in) :: name
    if (global_all(pass)) then
      if (is_IOP) write(output_unit,'(a)') 'Passed: ' //  name
    else
      status = 1
      if (is_IOP) write(output_unit,'(a)') 'FAILED: ' //  name
    end if
  end subroutine

  ! Distribute a scalar value to each process
  subroutine dist_scalar
    integer, allocatable :: asrc(:)
    integer :: n, adest
    asrc = [(n, n=1,npe)]
    adest = this_pe
    block
      integer(int32) :: src(npe), dest
      src = asrc
      dest = 0
      call distribute(src, dest)
      pass = (dest == adest)
      call write_result(pass, 'dist_scalar_int32')
    end block
    block
      integer(int64) :: src(npe), dest
      src = asrc
      dest = 0
      call distribute(src, dest)
      pass = (dest == adest)
      call write_result(pass, 'dist_scalar_int64')
    end block
    block
      real(real32) :: src(npe), dest
      src = asrc
      dest = 0
      call distribute(src, dest)
      pass = (dest == adest)
      call write_result(pass, 'dist_scalar_real32')
    end block
    block
      real(real64) :: src(npe), dest
      src = asrc
      dest = 0
      call distribute(src, dest)
      pass = (dest == adest)
      call write_result(pass, 'dist_scalar_real64')
    end block
  end subroutine

  ! generic rank-1 array case
  subroutine dist_rank1
    integer, allocatable :: asrc(:), adest(:)
    integer :: j, n
    if (is_IOP) then
      asrc = [(j, j=1, (npe*(npe+1))/2)]
    else
      allocate(asrc(0))
    end if
    n = (this_pe*(this_pe-1))/2
    adest = [(n+j, j=1,this_pe)]
    block
      integer(int8), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank1_int8')
    end block
    block
      integer(int32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank1_int32')
    end block
    block
      integer(int64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank1_int64')
    end block
    block
      real(real32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank1_real32')
    end block
    block
      real(real64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank1_real64')
    end block
  end subroutine

  ! Rank-1 array case with a 0-sized vector
  subroutine dist_rank1_zero
    integer, allocatable :: asrc(:), adest(:)
    integer :: j, n
    if (is_IOP) then
      asrc = [(j, j=1, (npe*(npe-1))/2)]
    else
      allocate(asrc(0))
    end if
    n = ((this_pe-1)*(this_pe-2))/2
    adest = [(n+j, j=1,this_pe-1)]
    block
      integer(int8), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank1_zero_int8')
    end block
    block
      integer(int32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank1_zero_int32')
    end block
    block
      integer(int64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank1_zero_int64')
    end block
    block
      real(real32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank1_zero_real32')
    end block
    block
      real(real64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank1_zero_real64')
    end block
  end subroutine

  ! generic rank-2 array case
  subroutine dist_rank2
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j, n
    if (is_IOP) then
      allocate(asrc(2,(npe*(npe+1))/2))
      asrc(1,:) = [(j, j=1, (npe*(npe+1))/2)]
      asrc(2,:) = -asrc(1,:)
    else
      allocate(asrc(2,0))
    end if
    n = (this_pe*(this_pe-1))/2
    allocate(adest(2,this_pe))
    adest(1,:) = [(n+j, j=1,this_pe)]
    adest(2,:) = -adest(1,:)
    block
      integer(int8), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank2_int8')
    end block
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank2_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank2_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank2_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank2_real64')
    end block
  end subroutine

  ! Rank-2 array case with a 0-sized vector
  subroutine dist_rank2_zero
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j, n
    if (is_IOP) then
      n = (npe*(npe-1))/2
      allocate(asrc(2,n))
      asrc(1,:) = [(j, j=1,n)]
      asrc(2,:) = -asrc(1,:)
    else
      allocate(asrc(2,0))
    end if
    n = ((this_pe-1)*(this_pe-2))/2
    allocate(adest(2,this_pe-1))
    adest(1,:) = [(n+j, j=1,this_pe-1)]
    adest(2,:) = -adest(1,:)
    block
      integer(int8), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank2_zero_int8')
    end block
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank2_zero_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank2_zero_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank2_zero_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank2_zero_real64')
    end block
  end subroutine

  ! rank-2 array section case
  subroutine dist_array_section
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j, n
    if (is_IOP) then
      n = (npe*(npe+1))/2
      allocate(asrc(3,2*n), source=-1)
      asrc(1,1::2) = [(j, j=1,n)]
      asrc(3,1::2) = -asrc(1,1::2)
    else
      allocate(asrc(3,0))
    end if
    n = (this_pe*(this_pe-1))/2
    allocate(adest(3,2*this_pe), source=0)
    adest(1,1::2) = [(n+j, j=1,this_pe)]
    adest(3,1::2) = -adest(1,1::2)
    block
      integer(int8), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0_int8)
      call distribute(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == adest)
      call write_result(pass, 'dist_array_section_int8')
    end block
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0_int32)
      call distribute(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == adest)
      call write_result(pass, 'dist_array_section_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0_int64)
      call distribute(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == adest)
      call write_result(pass, 'dist_array_section_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0.0_real32)
      call distribute(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == adest)
      call write_result(pass, 'dist_array_section_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0.0_real64)
      call distribute(src(1::2,1::2), dest(1::2,1::2))
      pass = all(dest == adest)
      call write_result(pass, 'dist_array_section_real64')
    end block
  end subroutine

  ! generic rank-3 array case
  subroutine dist_rank3
    integer, allocatable :: asrc(:,:,:), adest(:,:,:)
    integer :: j, n
    if (is_IOP) then
      allocate(asrc(2,2,(npe*(npe+1))/2))
      asrc(1,1,:) = [(j, j=1, (npe*(npe+1))/2)]
      asrc(2,1,:) = -asrc(1,1,:)
      asrc(1,2,:) = 2*asrc(1,1,:)
      asrc(2,2,:) = 2*asrc(2,1,:)
    else
      allocate(asrc(2,2,0))
    end if
    n = (this_pe*(this_pe-1))/2
    allocate(adest(2,2,this_pe))
    adest(1,1,:) = [(n+j, j=1,this_pe)]
    adest(2,1,:) = -adest(1,1,:)
    adest(1,2,:) = 2*adest(1,1,:)
    adest(2,2,:) = 2*adest(2,1,:)
    block
      integer(int8), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank3_int8')
    end block
    block
      integer(int32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank3_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank3_int64')
    end block
    block
      real(real32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank3_real32')
    end block
    block
      real(real64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank3_real64')
    end block
  end subroutine

  ! Rank-3 array case with a 0-sized vector
  subroutine dist_rank3_zero
    integer, allocatable :: asrc(:,:,:), adest(:,:,:)
    integer :: j, n
    if (is_IOP) then
      n = (npe*(npe-1))/2
      allocate(asrc(2,2,n))
      asrc(1,1,:) = [(j, j=1,n)]
      asrc(2,1,:) = -asrc(1,1,:)
      asrc(1,2,:) = 2*asrc(1,1,:)
      asrc(2,2,:) = 2*asrc(2,1,:)
    else
      allocate(asrc(2,2,0))
    end if
    n = ((this_pe-1)*(this_pe-2))/2
    allocate(adest(2,2,this_pe-1))
    adest(1,1,:) = [(n+j, j=1,this_pe-1)]
    adest(2,1,:) = -adest(1,1,:)
    adest(1,2,:) = 2*adest(1,1,:)
    adest(2,2,:) = 2*adest(2,1,:)
    block
      integer(int8), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank3_zero_int8')
    end block
    block
      integer(int32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank3_zero_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank3_zero_int64')
    end block
    block
      real(real32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank3_zero_real32')
    end block
    block
      real(real64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call distribute(src, dest)
      pass = all(dest == adest)
      call write_result(pass, 'dist_rank3_zero_real64')
    end block
  end subroutine

  subroutine dist_log_scalar
    logical, allocatable :: src(:)
    logical :: dest
    if (is_IOP) then
      allocate(src(npe), source=.true.)
    else
      allocate(src(0))
    end if
    dest = .false.
    call distribute(src, dest)
    pass = dest
    call write_result(pass, 'dist_log_scalar')
  end subroutine

  subroutine dist_log_rank1
    logical, allocatable :: src(:), dest(:)
    integer :: n
    if (is_IOP) then
      n = (npe*(npe+1))/2
      allocate(src(n), source=.true.)
    else
      allocate(src(0))
    end if
    allocate(dest(this_pe), source=.false.)
    call distribute(src, dest)
    pass = all(dest)
    call write_result(pass, 'dist_log_rank1')
  end subroutine

  subroutine dist_log_rank2
    logical, allocatable :: src(:,:), dest(:,:)
    integer :: n
    if (is_IOP) then
      n = (npe*(npe+1))/2
      allocate(src(2,n))
      src(1,:) = .true.
      src(2,:) = .false.
    else
      allocate(src(2,0))
    end if
    allocate(dest(2,this_pe))
    dest(1,:) = .false.
    dest(2,:) = .true.
    call distribute(src, dest)
    pass = all(dest(1,:) .and. .not.dest(2,:))
    call write_result(pass, 'dist_log_rank2')
  end subroutine

  subroutine dist_log_rank3
    logical, allocatable :: src(:,:,:), dest(:,:,:)
    integer :: n
    if (is_IOP) then
      n = (npe*(npe+1))/2
      allocate(src(2,2,n))
      src(1,1,:) = .true.
      src(2,1,:) = .false.
      src(1,2,:) = .false.
      src(2,2,:) = .true.
    else
      allocate(src(2,2,0))
    end if
    allocate(dest(2,2,this_pe))
    dest(1,1,:) = .false.
    dest(2,1,:) = .true.
    dest(1,2,:) = .true.
    dest(2,2,:) = .false.
    call distribute(src, dest)
    pass = all(dest(1,1,:) .and. .not.dest(2,1,:) .and. .not.dest(1,2,:) .and. dest(2,2,:))
    call write_result(pass, 'dist_log_rank3')
  end subroutine

end program