!! Unit Tests for INDEX_MAP Off-Process Gather Procedures
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

  use,intrinsic :: iso_fortran_env, only: int32, real32, real64, output_unit
  use index_map_type
#ifndef USE_CAF
  use mpi
#endif
  implicit none

  logical :: is_root
  integer :: my_rank, nproc, bsize, ierr, status
  integer, allocatable :: offp_index(:)
  type(index_map) :: imap

#ifdef USE_CAF
  my_rank = this_image()-1
  nproc = num_images()
#else
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
#endif
  is_root = (my_rank == 0)

  if (nproc /= 4) then
#ifndef USE_CAF
    call MPI_Finalize(ierr)
#endif
    if (is_root) write(output_unit,'(a)') 'Test must be run using 4 processes'
    error stop 1
  end if
  if (is_root) write(output_unit,'(a,i0,a)') 'Using ', nproc, ' processes'

  bsize = 3
  select case (my_rank)
  case (0)
    offp_index = [4,7]
  case (1)
    offp_index = [7,10,11]
  case (2)
    offp_index = [4,5,11,12]
  case (3)
    offp_index = [integer::]
  end select

  call imap%init(bsize, offP_index)

  status = 0
  call test_imap
  call test_rank1
  call test_rank2
  call test_rank3
  call test_log

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
      !sync all
      call co_min(n)
      !sync all
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

  subroutine test_imap
    integer :: j
    logical :: pass
    select case (my_rank)
    case (0)
      pass = (imap%onp_size==3) .and. (imap%local_size==5)
      pass = pass .and. all(imap%global_index([(j,j=1,5)]) == [1,2,3,4,7])
      call write_result(pass, 'test_imap')
    case (1)
      pass = (imap%onp_size==3) .and. (imap%local_size==6)
      pass = pass .and. all(imap%global_index([(j,j=1,6)]) == [4,5,6,7,10,11])
      call write_result(pass, 'test_imap')
    case (2)
      pass = (imap%onp_size==3) .and. (imap%local_size==7)
      pass = pass .and. all(imap%global_index([(j,j=1,7)]) == [7,8,9,4,5,11,12])
      call write_result(pass, 'test_imap')
    case (3)
      pass = (imap%onp_size==3) .and. (imap%local_size==3)
      pass = pass .and. all(imap%global_index([(j,j=1,3)]) == [10,11,12])
      call write_result(pass, 'test_imap')
    end select
  end subroutine

  ! We pad the array with an extra element to verify that its value is
  ! preserved by the gather operation. The wrong intent on the on-process
  ! array will result in it being overwritten when using NAG's -nan option.

  subroutine test_rank1
    integer ::j, input(imap%local_size+1), output(imap%local_size+1)
    input = 99
    do j = 1, imap%onp_size
      input(j) = imap%global_index(j)
    end do
    output = 99
    do j = 1, imap%local_size
      output(j) = imap%global_index(j)
    end do
    block
      integer(int32), allocatable :: array(:)
      array = input
      call imap%gather_offp(array)
      call write_result(all(array == output), 'test_rank1_int32')
    end block
    block
      real(real32), allocatable :: array(:)
      array = input
      call imap%gather_offp(array)
      call write_result(all(array == output), 'test_rank1_real32')
    end block
    block
      real(real64), allocatable :: array(:)
      array = input
      call imap%gather_offp(array)
      call write_result(all(array == output), 'test_rank1_real64')
    end block
    block
      complex(real32), allocatable :: array(:)
      array = (1,-1)*input
      call imap%gather_offp(array)
      call write_result(all(array == (1,-1)*output), 'test_rank1_complex32')
    end block
    block
      complex(real64), allocatable :: array(:)
      array = (1,-1)*input
      call imap%gather_offp(array)
      call write_result(all(array == (1,-1)*output), 'test_rank1_complex64')
    end block
  end subroutine

  subroutine test_rank2
    integer :: j, input(2,imap%local_size+1), output(2,imap%local_size+1)
    input = 99
    do j = 1, imap%onp_size
      input(1,j) = imap%global_index(j)
    end do
    input(2,:) = -input(1,:)
    output = 99
    do j = 1, imap%local_size
      output(1,j) = imap%global_index(j)
    end do
    output(2,:) = -output(1,:)
    block
      integer(int32), allocatable :: array(:,:)
      array = input
      call imap%gather_offp(array)
      call write_result(all(array == output), 'test_rank2_int32')
    end block
    block
      real(real32), allocatable :: array(:,:)
      array = input
      call imap%gather_offp(array)
      call write_result(all(array == output), 'test_rank2_real32')
    end block
    block
      real(real64), allocatable :: array(:,:)
      array = input
      call imap%gather_offp(array)
      call write_result(all(array == output), 'test_rank2_real64')
    end block
    block
      complex(real32), allocatable :: array(:,:)
      array = (1,-1)*input
      call imap%gather_offp(array)
      call write_result(all(array == (1,-1)*output), 'test_rank2_complex32')
    end block
    block
      complex(real64), allocatable :: array(:,:)
      array = (1,-1)*input
      call imap%gather_offp(array)
      call write_result(all(array == (1,-1)*output), 'test_rank2_complex64')
    end block
  end subroutine

  subroutine test_rank3
    integer :: j, input(2,2,imap%local_size+1), output(2,2,imap%local_size+1)
    input = 99
    do j = 1, imap%onp_size
      input(1,1,j) = imap%global_index(j)
    end do
    input(2,1,:) = -input(1,1,:)
    input(:,2,:) = 2*input(:,1,:)
    output = 99
    do j = 1, imap%local_size
      output(1,1,j) = imap%global_index(j)
    end do
    output(2,1,:) = -output(1,1,:)
    output(:,2,:) = 2*output(:,1,:)
    block
      integer(int32), allocatable :: array(:,:,:)
      array = input
      call imap%gather_offp(array)
      call write_result(all(array == output), 'test_rank3_int32')
    end block
    block
      real(real32), allocatable :: array(:,:,:)
      array = input
      call imap%gather_offp(array)
      call write_result(all(array == output), 'test_rank3_real32')
    end block
    block
      real(real64), allocatable :: array(:,:,:)
      array = input
      call imap%gather_offp(array)
      call write_result(all(array == output), 'test_rank3_real64')
    end block
    block
      complex(real32), allocatable :: array(:,:,:)
      array = (1,-1)*input
      call imap%gather_offp(array)
      call write_result(all(array == (1,-1)*output), 'test_rank3_complex32')
    end block
    block
      complex(real64), allocatable :: array(:,:,:)
      array = (1,-1)*input
      call imap%gather_offp(array)
      call write_result(all(array == (1,-1)*output), 'test_rank3_complex64')
    end block
  end subroutine

  subroutine test_log
    integer :: j
    logical :: input(imap%local_size+1), output(imap%local_size+1)
    input = .true.
    do j = 1, imap%local_size
      input(j) = (modulo(imap%global_index(j),2) == 0)
    end do
    do j = imap%onp_size+1, imap%local_size
      input(j) = .not.input(j)
    end do
    output = .true.
    do j = 1, imap%local_size
      output(j) = (modulo(imap%global_index(j),2) == 0)
    end do
    block
      logical, allocatable :: array(:)
      array = input
      call imap%gather_offp(array)
      call write_result(all(array .eqv. output), 'test_log_rank1')
    end block
    block
      logical :: array(2,imap%local_size+1)
      array(1,:) = input
      array(2,:) = .not.input
      call imap%gather_offp(array)
      call write_result(all((array(1,:) .eqv. output) .and. (array(2,:) .neqv. output)), 'test_log_rank2')
    end block
    block
      logical :: array(2,2,imap%local_size+1)
      array(1,1,:) = input
      array(2,1,:) = .not.input
      array(:,2,:) = .not.array(:,1,:)
      call imap%gather_offp(array)
      call write_result( &
          all((array(1,1,:) .eqv. output) .and. (array(2,1,:) .neqv. output) .and. &
              (array(1,2,:) .neqv. output) .and. (array(2,2,:) .eqv. output)), 'test_log_rank3')
    end block
  end subroutine

end program
