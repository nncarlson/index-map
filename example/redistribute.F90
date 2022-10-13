!! REDISTRIBUTE
!!
!! This example illustrates how one can use the INDEX_MAP derived type to
!! redistribute an array that is distributed according to one partitioning
!! of its index set to another array that is distributed according to another
!! partitioning of the index set.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

program redistribute

#ifndef USE_CAF
  use mpi
#endif
  use index_map_type
  implicit none

  integer :: ierr, nproc, my_rank, bsize, j, n, stat
  real, allocatable :: a1(:), a2(:), a1_global(:), a2_global(:)
  type(index_map) :: imap1, imap2

#ifdef USE_CAF
  nproc = num_images()
  my_rank = this_image() - 1
  if (my_rank == 0) write(*,'(a,i0,a)') 'Running with ', nproc, ' CAF images'
#else
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
  if (my_rank == 0) write(*,'(a,i0,a)') 'Running with ', nproc, ' MPI ranks'
#endif

  !! Start with a distributed array A1 with respect to IMAP1.
  bsize = 2*(my_rank+1)
  call imap1%init(bsize)
  allocate(a1(imap1%onp_size))
  do j = 1, size(a1)
    a1(j) = imap1%global_index(j)
  end do

  !! Consider a different distribution IMAP2 of the same global index set.
  bsize = 2*(nproc - my_rank)
  call imap2%init(bsize)
  ASSERT(imap2%global_size == imap1%global_size)

  !! The object now is to redistribute the distributed values of A1
  !! to the array A2 that is distributed with respect to IMAP2.
  allocate(a2(imap2%onp_size))
  block
    type(index_map) :: imap1_tmp
    integer, allocatable :: p(:)
    real, allocatable :: a1_offp(:)

    !! We need a copy of IMAP1 that we can add off-process indices to.
    call imap1_tmp%init(imap1%onp_size)

    !! To start, P maps local indices with respect to IMAP2 to global indices. 
    allocate(p(imap2%onp_size))
    do j = 1, imap2%onp_size
      p(j) = imap2%global_index(j)
    end do
    
    !! Remap the global P index values to local indices with respect to IMAP1.
    call imap1_tmp%localize_index_array(p)
    !! P now maps IMAP2 local indices to IMAP1 local indices some of which are
    !! are off-process that were added by the call to LOCALIZE_INDEX_ARRAY.

    !! Gather the off-process values of A1.
    allocate(a1_offp(imap1_tmp%onp_size+1:imap1_tmp%local_size))
    call imap1_tmp%gather_offp(a1, a1_offp)

    !! Finally we can fill the distributed A2 array with purely local copies.
    do j = 1, size(a2)
      n = p(j)
      if (n > imap1_tmp%onp_size) then
        a2(j) = a1_offp(n)
      else
        a2(j) = a1(n)
      end if
    end do
  end block

  !! Verify that the result is correct
  n = merge(imap1%global_size, 0, my_rank==0)
  allocate(a1_global(n), a2_global(n))
  call imap1%collate(a1, a1_global)
  call imap2%collate(a2, a2_global)
  if (my_rank == 0) stat = merge(0, 1, all(a1_global == a2_global))
#ifdef USE_CAF
  call co_broadcast(stat, source_image=1)
#else
  call MPI_Bcast(stat, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Finalize(ierr)
#endif
  if (stat == 0) then
    if (my_rank == 0) write(*,'(a)') 'Success!'
  else
    if (my_rank == 0) write(*,'(a)') 'Failure!'
    error stop 1
  end if

end program
