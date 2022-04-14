!! Unit Tests for INDEX_MAP Distribute Procedures
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

#ifndef USE_CAF
  use mpi
#endif
  use index_map_type
  use,intrinsic :: iso_fortran_env
  implicit none

  integer :: ierr, my_rank, nproc, status
  logical :: is_root

#ifdef USE_CAF
  my_rank = this_image()-1
  nproc = num_images()
#else
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
#endif
  is_root = (my_rank == 0)
  if (is_root) write(output_unit,'(a,i0,a)') 'Using ', nproc, ' processes'

  status = 0
  call dist_rank1
  call dist_rank1_zero
  call dist_rank2
  call dist_rank2_zero
  call dist_array_section
  call dist_rank3
  call dist_rank3_zero
  call dist_log_rank1
  call dist_log_rank2
  call dist_log_rank3

#ifndef USE_CAF
  call MPI_Finalize(ierr)
#endif
  if (status /= 0) error stop 1

contains

  subroutine write_result(pass, name)
    logical, value :: pass
    character(*), intent(in) :: name
#ifdef USE_CAF
    block
      integer :: n
      n = merge(1, 0, pass)
      call co_min(n)
      pass = (n == 1)
    end block
#else
    call MPI_Allreduce(MPI_IN_PLACE, pass, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
#endif
    if (pass) then
      if (is_root) write(output_unit,'(a)') 'Passed: ' //  name
    else
      status = 1
      if (is_root) write(output_unit,'(a)') 'FAILED: ' //  name
    end if
  end subroutine

  ! generic rank-1 array case
  subroutine dist_rank1
    type(index_map) :: imap
    integer, allocatable :: asrc(:), adest(:)
    integer :: j
    call imap%init(1+my_rank)
    if (is_root) then
      asrc = [(j, j=1, imap%global_size)]
    else
      allocate(asrc(0))
    end if
    adest = [(j, j=imap%first_gid,imap%last_gid)]
    block
      integer(int32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank1_int32')
    end block
    block
      integer(int64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank1_int64')
    end block
    block
      real(real32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank1_real32')
    end block
    block
      real(real64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank1_real64')
    end block
  end subroutine

  ! Rank-1 array case with a 0-sized vector
  subroutine dist_rank1_zero
    type(index_map) :: imap
    integer, allocatable :: asrc(:), adest(:)
    integer :: j
    call imap%init(modulo(my_rank+2,nproc))
    if (is_root) then
      asrc = [(j, j=1, imap%global_size)]
    else
      allocate(asrc(0))
    end if
    adest = [(j, j=imap%first_gid,imap%last_gid)]
    block
      integer(int32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank1_zero_int32')
    end block
    block
      integer(int64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank1_zero_int64')
    end block
    block
      real(real32), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank1_zero_real32')
    end block
    block
      real(real64), allocatable :: src(:), dest(:)
      src = asrc
      allocate(dest(size(adest)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank1_zero_real64')
    end block
  end subroutine

  ! generic rank-2 array case
  subroutine dist_rank2
    type(index_map) :: imap
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(asrc(2,imap%global_size))
      asrc(1,:) = [(j, j=1, imap%global_size)]
      asrc(2,:) = -asrc(1,:)
    else
      allocate(asrc(2,0))
    end if
    allocate(adest(2,imap%onp_size))
    adest(1,:) = [(j, j=imap%first_gid,imap%last_gid)]
    adest(2,:) = -adest(1,:)
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank2_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank2_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank2_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank2_real64')
    end block
  end subroutine

  ! Rank-2 array case with a 0-sized vector
  subroutine dist_rank2_zero
    type(index_map) :: imap
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j
    call imap%init(modulo(my_rank+nproc-1,nproc))
    if (is_root) then
      allocate(asrc(2,imap%global_size))
      asrc(1,:) = [(j, j=1, imap%global_size)]
      asrc(2,:) = -asrc(1,:)
    else
      allocate(asrc(2,0))
    end if
    allocate(adest(2,imap%onp_size))
    adest(1,:) = [(j, j=imap%first_gid,imap%last_gid)]
    adest(2,:) = -adest(1,:)
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank2_zero_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank2_zero_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank2_zero_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(2,size(adest,2)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank2_zero_real64')
    end block
  end subroutine

  ! rank-2 array section case
  subroutine dist_array_section
    type(index_map) :: imap
    integer, allocatable :: asrc(:,:), adest(:,:)
    integer :: j
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(asrc(3,2*imap%global_size), source=127)
      asrc(1,1::2) = [(j, j=1, imap%global_size)]
      asrc(3,1::2) = -asrc(1,1::2)
    else
      allocate(asrc(3,0))
    end if
    allocate(adest(3,2*imap%onp_size), source=0)
    adest(1,1::2) = [(j, j=imap%first_gid,imap%last_gid)]
    adest(3,1::2) = -adest(1,1::2)
    block
      integer(int32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0_int32)
      call imap%distribute(src(1::2,1::2), dest(1::2,1::2))
      call write_result(all(dest == adest), 'dist_array_section_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0_int64)
      call imap%distribute(src(1::2,1::2), dest(1::2,1::2))
      call write_result(all(dest == adest), 'dist_array_section_int64')
    end block
    block
      real(real32), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0.0_real32)
      call imap%distribute(src(1::2,1::2), dest(1::2,1::2))
      call write_result(all(dest == adest), 'dist_array_section_real32')
    end block
    block
      real(real64), allocatable :: src(:,:), dest(:,:)
      src = asrc
      allocate(dest(3,size(adest,2)), source=0.0_real64)
      call imap%distribute(src(1::2,1::2), dest(1::2,1::2))
      call write_result(all(dest == adest), 'dist_array_section_real64')
    end block
  end subroutine

  ! generic rank-3 array case
  subroutine dist_rank3
    type(index_map) :: imap
    integer, allocatable :: asrc(:,:,:), adest(:,:,:)
    integer :: j
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(asrc(2,2,imap%global_size))
      asrc(1,1,:) = [(j, j=1, imap%global_size)]
      asrc(2,1,:) = -asrc(1,1,:)
      asrc(:,2,:) = 2*asrc(:,1,:)
    else
      allocate(asrc(2,2,0))
    end if
    allocate(adest(2,2,imap%onp_size))
    adest(1,1,:) = [(j, j=imap%first_gid,imap%last_gid)]
    adest(2,1,:) = -adest(1,1,:)
    adest(:,2,:) = 2*adest(:,1,:)
    block
      integer(int32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank3_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank3_int64')
    end block
    block
      real(real32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank3_real32')
    end block
    block
      real(real64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank3_real64')
    end block
  end subroutine

  ! Rank-3 array case with a 0-sized vector
  subroutine dist_rank3_zero
    type(index_map) :: imap
    integer, allocatable :: asrc(:,:,:), adest(:,:,:)
    integer :: j
    call imap%init(modulo(my_rank+nproc-1,nproc))
    if (is_root) then
      allocate(asrc(2,2,imap%global_size))
      asrc(1,1,:) = [(j, j=1, imap%global_size)]
      asrc(2,1,:) = -asrc(1,1,:)
      asrc(:,2,:) = 2*asrc(:,1,:)
    else
      allocate(asrc(2,2,0))
    end if
    allocate(adest(2,2,imap%onp_size))
    adest(1,1,:) = [(j, j=imap%first_gid,imap%last_gid)]
    adest(2,1,:) = -adest(1,1,:)
    adest(:,2,:) = 2*adest(:,1,:)
    block
      integer(int32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank3_zero_int32')
    end block
    block
      integer(int64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank3_zero_int64')
    end block
    block
      real(real32), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank3_zero_real32')
    end block
    block
      real(real64), allocatable :: src(:,:,:), dest(:,:,:)
      src = asrc
      allocate(dest(2,2,size(adest,3)))
      call imap%distribute(src, dest)
      call write_result(all(dest == adest), 'dist_rank3_zero_real64')
    end block
  end subroutine

  subroutine dist_log_rank1
    type(index_map) :: imap
    logical, allocatable :: src(:), dest(:)
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(src(imap%global_size), source=.true.)
    else
      allocate(src(0))
    end if
    allocate(dest(imap%onp_size), source=.false.)
    call imap%distribute(src, dest)
    call write_result(all(dest), 'dist_log_rank1')
  end subroutine

  subroutine dist_log_rank2
    type(index_map) :: imap
    logical, allocatable :: src(:,:), dest(:,:)
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(src(2,imap%global_size))
      src(1,:) = .true.
      src(2,:) = .false.
    else
      allocate(src(2,0))
    end if
    allocate(dest(2,imap%onp_size))
    dest(1,:) = .false.
    dest(2,:) = .true.
    call imap%distribute(src, dest)
    call write_result(all(dest(1,:) .and. .not.dest(2,:)), 'dist_log_rank2')
  end subroutine

  subroutine dist_log_rank3
    type(index_map) :: imap
    logical, allocatable :: src(:,:,:), dest(:,:,:)
    call imap%init(1+my_rank)
    if (is_root) then
      allocate(src(2,2,imap%global_size))
      src(1,1,:) = .true.
      src(2,1,:) = .false.
      src(:,2,:) = .not.src(:,1,:)
    else
      allocate(src(2,2,0))
    end if
    allocate(dest(2,2,imap%onp_size))
    dest(1,1,:) = .false.
    dest(2,1,:) = .true.
    dest(:,2,:) = .not.dest(:,1,:)
    call imap%distribute(src, dest)
    call write_result(all(dest(1,1,:) .and. .not.dest(2,1,:) .and. .not.dest(1,2,:) .and. dest(2,2,:)), 'dist_log_rank3')
  end subroutine

end program
