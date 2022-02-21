!! Unit Tests for INDEX_MAP Localize Procedures
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

  use,intrinsic :: iso_fortran_env, only: output_unit
  use index_map_type
  use mpi
  implicit none

  logical :: is_root
  integer :: my_rank, nproc, ierr, status

  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  is_root = (my_rank == 0)

  status = 0
  call test_basic
  call test_rank1
  call test_rank2
  call test_rank1_offp
  call test_rank1_domain_offp
  call test_struct

  call MPI_Finalize(ierr)
  if (status /= 0) error stop 1

contains

  subroutine write_result(pass, name)
    logical, value :: pass
    character(*), intent(in) :: name
    call MPI_Allreduce(MPI_IN_PLACE, pass, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    if (pass) then
      if (is_root) write(output_unit,'(a)') 'Passed: ' //  name
    else
      status = 1
      if (is_root) write(output_unit,'(a)') 'FAILED: ' //  name
    end if
  end subroutine

  ! No off-process references

  subroutine test_basic

    type(index_map) :: domain, range
    integer, allocatable :: g_index(:), l_index(:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:), l_output(:)
    integer :: j

    call domain%init(MPI_COMM_WORLD, 2)
    call range%init(MPI_COMM_WORLD, 4)

    g_index = [(2*j, j=1,merge(domain%global_size,0,is_root))]
    call domain%localize_index_array(g_index, range, l_index)
    call write_result(.not.allocated(range%offp_index), 'test_basic_no-offp')

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%distribute(g_input, l_input)

    g_output = g_input(g_index)
    allocate(l_output(domain%local_size))
    call domain%distribute(g_output, l_output)

    call write_result(all(l_output == l_input(l_index)), 'test_basic')

  end subroutine

  ! Rank-1 indexing array with off-process references
  ! No off-process in domain or range.
  ! Distributed index_map initialization.

  subroutine test_rank1

    type(index_map) :: domain, range
    integer, allocatable :: g_index(:), l_index(:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:), l_output(:)
    integer :: j

    call domain%init(MPI_COMM_WORLD, 1+my_rank)
    call range%init(MPI_COMM_WORLD, 2*(nproc-my_rank))

    g_index = [(2*j, j=1,merge(domain%global_size,0,is_root))]
    call domain%localize_index_array(g_index, range, l_index)

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%distribute(g_input, l_input)
    call range%gather_offp(l_input)

    g_output = g_input(g_index)
    allocate(l_output(domain%local_size))
    call domain%distribute(g_output, l_output)

    call write_result(all(l_output == l_input(l_index)), 'test_rank1')

  end subroutine

  ! Rank-2 indexing array with off-process references
  ! No off-process in domain or range.
  ! Distributed index_map initialization.

  subroutine test_rank2

    type(index_map) :: domain, range
    integer, allocatable :: g_index(:,:), l_index(:,:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:,:), l_output(:,:)
    integer :: j

    call domain%init(MPI_COMM_WORLD, 1+my_rank)
    call range%init(MPI_COMM_WORLD, 2*(nproc-my_rank))

    allocate(g_index(2,merge(domain%global_size,0,is_root)))
    g_index(1,:) = [(2*j, j=1,merge(domain%global_size,0,is_root))]
    g_index(2,:) = [(2*j-1, j=1,merge(domain%global_size,0,is_root))]
    call domain%localize_index_array(g_index, range, l_index)

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%distribute(g_input, l_input)
    call range%gather_offp(l_input)

    allocate(g_output, mold=g_index)
    g_output(1,:) = g_input(g_index(1,:))
    g_output(2,:) = g_input(g_index(2,:))
    allocate(l_output, mold=l_index)
    call domain%distribute(g_output, l_output)

    call write_result(all(l_output(1,:) == l_input(l_index(1,:)) .and. &
                          l_output(2,:) == l_input(l_index(2,:))), 'test_rank2')

  end subroutine

  ! Rank-1 indexing array with off-process references
  ! No off-process in domain, but range has pre-existing off-process
  ! Exercises distributed index_map initialization

  subroutine test_rank1_offp

    type(index_map) :: domain, range
    integer, allocatable :: g_index(:), l_index(:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:), l_output(:)
    integer, allocatable :: offp_index(:)
    integer :: j, bsize

    bsize = 3
    offp_index = [1+modulo((my_rank+1)*bsize, nproc*bsize)]
    call range%init(MPI_COMM_WORLD, bsize, offp_index)
    call domain%init(MPI_COMM_WORLD, 1)

    g_index = [(1+modulo(j,nproc)*bsize, j=1,nproc)]
    call domain%localize_index_array(g_index, range, l_index)

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%distribute(g_input, l_input)
    call range%gather_offp(l_input)

    g_output = g_input(g_index)
    allocate(l_output(domain%local_size))
    call domain%distribute(g_output, l_output)

    call write_result(all(l_output == l_input(l_index)), 'test_rank1_offp')

  end subroutine

  ! Rank-1 indexing array with off-process references
  ! Domain has off-process elements; no off-process in range.
  ! Distributed and root index_map initialization.

  subroutine test_rank1_domain_offp

    type(index_map) :: domain, range
    integer, allocatable :: g_index(:), l_index(:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:), l_output(:)
    integer, allocatable :: bsizes(:), offp_counts(:), offp_indices(:)
    integer :: j, n

    if (is_root) then
      allocate(bsizes(nproc), offp_counts(nproc), offp_indices(nproc))
      bsizes = [(j, j=1,nproc)]
      offp_counts = 1
      n = (nproc*(nproc+1))/2
      offp_indices = [(1+modulo((j*(j+1))/2,n), j=1,nproc)]
    else
      allocate(bsizes(0), offp_counts(0), offp_indices(0))
    end if

    call domain%init(MPI_COMM_WORLD, bsizes, offp_counts, offp_indices)
    call range%init(MPI_COMM_WORLD, 2*(nproc-my_rank))

    g_index = [(2*j, j=1,merge(domain%global_size,0,is_root))]
    call domain%localize_index_array(g_index, range, l_index)

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%distribute(g_input, l_input)
    call range%gather_offp(l_input)

    g_output = g_input(g_index)
    allocate(l_output(domain%local_size))
    call domain%distribute(g_output, l_output)
    call domain%gather_offp(l_output)

    call write_result(all(l_output == l_input(l_index)), 'test_rank1_domain_offp')

  end subroutine

  ! Indexing array structure with off-process reference and
  ! domain with off-process

  subroutine test_struct

    type(index_map) :: domain, range
    integer, allocatable :: g_count(:), g_index(:), l_count(:), l_index(:)
    integer, allocatable :: g_input(:), l_input(:), g_output(:), l_output(:)

    integer :: j, n, m

    m = (nproc*(nproc+1))/2
    n = 1 + modulo(((my_rank+1)*(my_rank+2))/2, m)
    call domain%init(MPI_COMM_WORLD, my_rank+1, offp_index=[n])
    call range%init(MPI_COMM_WORLD, nproc-my_rank)

    g_count = [(modulo(j,3), j=1,merge(domain%global_size,0,is_root))]
    g_index = [(1+modulo(j,range%global_size), j=1,sum(g_count))]
    call domain%localize_index_array(g_count, g_index, range, l_count, l_index)

    g_input = [(j, j=1,merge(range%global_size,0,is_root))]
    allocate(l_input(range%local_size))
    call range%distribute(g_input, l_input)
    call range%gather_offp(l_input)

    g_output = g_input(g_index)
    allocate(l_output(domain%local_size))
    call domain%distribute(g_output, l_output)
    call domain%gather_offp(l_output)

    call write_result(all(l_output == l_input(l_index)), 'test_struct')

  end subroutine

end program
