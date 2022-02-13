!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2022  Neil N. Carlson
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

module index_map_type

  use mpi
  use,intrinsic :: iso_fortran_env, only: int32, real64
  implicit none
  private

  type, public :: index_map
    integer :: comm ! PROTECTED
    integer :: onp_size=0, offp_size=0, local_size=0  ! PROTECTED
    integer :: global_size=0, first=0, last=0 ! PROTECTED; int64 option?
    integer, allocatable :: offp_index(:)  ! int64 option?
    integer, allocatable, private :: offp_counts(:), offp_displs(:)
    integer, allocatable, private :: onp_index(:), onp_counts(:), onp_displs(:)
  contains
    generic   :: init => init1, init2, init3
    procedure :: add_offp_index
    generic   :: gather => gather1_int32_r1, gather2_int32_r1, &
                           gather1_int32_r2, gather2_int32_r2, &
                           gather1_int32_r3, gather2_int32_r3, &
                           gather1_real64_r1, gather2_real64_r1, &
                           gather1_real64_r2, gather2_real64_r2, &
                           gather1_real64_r3, gather2_real64_r3
    generic   :: scatter_sum => scatter1_sum_int32_r1, scatter2_sum_int32_r1
    generic   :: scatter_min => scatter1_min_int32_r1, scatter2_min_int32_r1
    generic   :: scatter_max => scatter1_max_int32_r1, scatter2_max_int32_r1
    generic   :: scatter_or  => scatter1_or_l_r1,  scatter2_or_l_r1
    generic   :: scatter_and => scatter1_and_l_r1, scatter2_and_l_r1
    procedure :: global_index
    procedure, private :: init1, init2, init3
    procedure, private ::  gather1_int32_r1, gather2_int32_r1, &
                           gather1_int32_r2, gather2_int32_r2, &
                           gather1_int32_r3, gather2_int32_r3, &
                           gather1_real64_r1, gather2_real64_r1, &
                           gather1_real64_r2, gather2_real64_r2, &
                           gather1_real64_r3, gather2_real64_r3
    procedure, private :: scatter1_sum_int32_r1, scatter2_sum_int32_r1
    procedure, private :: scatter1_min_int32_r1, scatter2_min_int32_r1
    procedure, private :: scatter1_max_int32_r1, scatter2_max_int32_r1
    procedure, private :: scatter1_or_l_r1,  scatter2_or_l_r1
    procedure, private :: scatter1_and_l_r1, scatter2_and_l_r1
  end type

  interface
    module subroutine gather1_int32_r1(this, local_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gather2_int32_r1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(in) :: onp_data(:)
      integer(int32), intent(out) :: offp_data(:)
    end subroutine
    module subroutine gather1_int32_r2(this, local_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gather2_int32_r2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(in) :: onp_data(:,:)
      integer(int32), intent(out) :: offp_data(:,:)
    end subroutine
    module subroutine gather1_int32_r3(this, local_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gather2_int32_r3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(in) :: onp_data(:,:,:)
      integer(int32), intent(out) :: offp_data(:,:,:)
    end subroutine
    module subroutine gather1_real64_r1(this, local_data)
      class(index_map), intent(in) :: this
      real(real64), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gather2_real64_r1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(real64), intent(in) :: onp_data(:)
      real(real64), intent(out) :: offp_data(:)
    end subroutine
    module subroutine gather1_real64_r2(this, local_data)
      class(index_map), intent(in) :: this
      real(real64), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gather2_real64_r2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(real64), intent(in) :: onp_data(:,:)
      real(real64), intent(out) :: offp_data(:,:)
    end subroutine
    module subroutine gather1_real64_r3(this, local_data)
      class(index_map), intent(in) :: this
      real(real64), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gather2_real64_r3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(real64), intent(in) :: onp_data(:,:,:)
      real(real64), intent(out) :: offp_data(:,:,:)
    end subroutine
  end interface

  interface
    module subroutine scatter1_sum_int32_r1(this, local_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scatter2_sum_int32_r1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(inout) :: onp_data(:)
      integer(int32), intent(in) :: offp_data(:)
    end subroutine
    module subroutine scatter1_min_int32_r1(this, local_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scatter2_min_int32_r1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(inout) :: onp_data(:)
      integer(int32), intent(in) :: offp_data(:)
    end subroutine
    module subroutine scatter1_max_int32_r1(this, local_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scatter2_max_int32_r1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(int32), intent(inout) :: onp_data(:)
      integer(int32), intent(in) :: offp_data(:)
    end subroutine
    module subroutine scatter1_or_l_r1(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:)
    end subroutine
    module subroutine scatter2_or_l_r1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: onp_data(:)
      logical, intent(in) :: offp_data(:)
    end subroutine
    module subroutine scatter1_and_l_r1(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:)
    end subroutine
    module subroutine scatter2_and_l_r1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: onp_data(:)
      logical, intent(in) :: offp_data(:)
    end subroutine
  end interface

!  interface localize_index_array
!    module subroutine localize_index_array1(g_index, domain, range, l_index, offp_index)
!      integer, intent(in) :: g_index(:)
!      class(index_map), intent(in) :: domain, range
!      integer, allocatable :: l_index(:), offP_index(:)
!    end subroutine
!  end interface
  
contains

  !! Each rank supplied its block size
  subroutine init1(this, comm, bsize)

    class(index_map), intent(out) :: this
    integer, intent(in) :: comm
    integer, intent(in) :: bsize

    integer :: np, ierr

    this%comm = comm

    call MPI_Comm_size(comm, np, ierr)

    this%onp_size = bsize
    this%offp_size = 0
    this%local_size = this%onp_size + this%offp_size
    call MPI_Scan(bsize, this%last, 1, MPI_INTEGER, MPI_SUM, this%comm, ierr)
    this%first = this%last - this%onp_size + 1
    this%global_size = this%last
    call MPI_Bcast(this%global_size, 1, MPI_INTEGER, np-1, this%comm, ierr)

  end subroutine init1

  !! One root rank has an array of rank block sizes
  subroutine init2(this, comm, bsizes, root)
    class(index_map), intent(out) :: this
    integer, intent(in) :: comm
    integer, intent(in) :: bsizes(:)
    integer, intent(in), optional :: root
    integer :: np, ierr, my_rank, root_, bsize
    call MPI_Comm_rank(comm, my_rank, ierr)
    root_ = 0
    if (present(root)) root_ = root
    if (my_rank == root_) then
      call MPI_Comm_size(comm, np, ierr)
      INSIST(size(bsizes) == np)
    end if
    call MPI_Scatter(bsizes, 1, MPI_INTEGER, bsize, 1, MPI_INTEGER, root_, comm, ierr)
    call init1(this, comm, bsize)
  end subroutine init2

  !! Each rank supplied with its block size and list of off-process indices
  subroutine init3(this, comm, bsize, offp_index)
    class(index_map), intent(out) :: this
    integer, intent(in) :: comm
    integer, intent(in) :: bsize, offp_index(:)
    call init1(this, comm, bsize)
    call add_offp_index(this, offp_index)
  end subroutine init3

  !! Add off-process indices to an already initialized index map
  subroutine add_offp_index(this, offp_index)

    class(index_map), intent(inout) :: this
    integer, intent(in) :: offp_index(:)

    integer :: ierr, np, my_rank, rank, i, j, j1, n
    integer, allocatable :: last(:), offp_rank(:)
    integer, allocatable :: onp_count(:), offp_count(:)
    integer, allocatable :: onp_ranks(:), offp_ranks(:)
    integer :: new_comm
    
    ASSERT(minval(offP_index) >= 1)
    ASSERT(maxval(offP_index) <= this%global_size)
    ASSERT(all((offp_index < this%first) .or. (offp_index > this%last)))

    !TODO? Add code to sort offp_index and remove duplicates
    ASSERT(all(offp_index(2:) > offp_index(:size(offp_index)-1)))

    !TODO? Allow extending an existing %offp_index
    ASSERT(.not.allocated(this%offp_index))

    this%offp_index = offp_index
    this%offp_size  = size(this%offp_index)
    this%local_size = this%onp_size + this%offp_size

    call MPI_Comm_size(this%comm, np, ierr)
    call MPI_Comm_rank(this%comm, my_rank, ierr)

    !! Determine which rank owns each off-process index (OFFP_RANK).
    !! OFFP_RANK will be ordered if OFFP_INDEX is ordered (we need this later).
    allocate(last(0:np-1), offp_rank(this%offp_size))
    call MPI_Allgather(this%last, 1, MPI_INTEGER, last, 1, MPI_INTEGER, this%comm, ierr)
    rank = 0
    do j = 1, size(this%offp_index)
      do while (this%offp_index(j) > last(rank))
        rank = rank + 1
        ASSERT(rank < np)
      end do
      ASSERT(rank /= my_rank)
      offp_rank(j) = rank
    end do
    deallocate(last)

    !! Count the number of our off-process indices that are owned by each
    !! rank (OFFP_COUNT). This relies on the OFFP_RANK array being ordered.
    allocate(offp_count(0:np-1))
    j1 = 1
    do rank = 0, np-1
      do j = j1, size(offp_rank)
        if (offp_rank(j) > rank) exit
      end do
      offp_count(rank) = j - j1
      j1 = j
    end do
    deallocate(offp_rank)

    !! Distribute the number of our off-process indices that are owned by each
    !! rank to that rank. The result is the number of our on-process indices
    !! that are off-process on each rank (ONP_COUNT).
    allocate(onp_count(0:np-1))
    call MPI_Alltoall(offp_count, 1, MPI_INTEGER, onp_count, 1, MPI_INTEGER, this%comm, ierr)

    !! We now know which ranks need to communicate with which ranks in order
    !! to exchange data between off-process indices and their corresponding
    !! on-process indices. This all that is needed to define a virtual graph
    !! process topology.

    !! Compress the OFFP_COUNT array and generate the list of ranks that own
    !! our off-process indices: OFFP_RANKS is the list of incoming neighbors.
    n = count(offp_count > 0)
    allocate(offp_ranks(n), this%offp_counts(n))
    n = 0
    do rank = 0, np-1
      if (offp_count(rank) > 0) then
        n = n + 1
        offp_ranks(n) = rank
        this%offp_counts(n) = offp_count(rank)
      end if
    end do
    deallocate(offp_count)

    !! Compress the ONP_COUNT array and generate the list of ranks with off-
    !! process indices owned by us: ONP_RANKS is the list of outgoing neighbors.
    n = count(onp_count > 0)
    allocate(onp_ranks(n), this%onp_counts(n))
    n = 0
    do rank = 0, np-1
      if (onp_count(rank) > 0) then
        n = n + 1
        onp_ranks(n) = rank
        this%onp_counts(n) = onp_count(rank)
      end if
    end do
    deallocate(onp_count)

    !! Create the virtual topology. Here we use the OFFP_COUNTS and ONP_COUNTS
    !! as the edge weights (used if rank reordering is enabled). An alternative
    !! is to replace them with MPI_UNWEIGHTED.  Rank reordering is disabled.
    call MPI_Dist_graph_create_adjacent(this%comm, &
        size(offp_ranks), offp_ranks, this%offp_counts, &
        size(onp_ranks),  onp_ranks,  this%onp_counts, &
        MPI_INFO_NULL, .false., new_comm, ierr)
    this%comm = new_comm
    deallocate(onp_ranks, offp_ranks)

    !! The components %OFFP_COUNTS, %OFFP_DISPLS, %ONP_COUNTS and %ONP_DIPSLS
    !! initialized here are meant for use with MPI_Neighbor_alltoallv. The
    !! component %ONP_INDEX will be used to fill the on-process buffer; the
    !! corresponding off-process buffer is off-process data array itself.

    !! Generate displacements into the on-process buffer for the start of the
    !! data for each neighbor rank.
    allocate(this%onp_displs, mold=this%onp_counts)
    if (size(this%onp_displs) > 0) then
      this%onp_displs(1) = 0
      do i = 2, size(this%onp_displs)
        this%onp_displs(i) = this%onp_displs(i-1) + this%onp_counts(i-1)
      end do
    end if

    !! Generate displacements into the off-process buffer for the start of
    !! the data for each neighbor rank.
    allocate(this%offp_displs, mold=this%offp_counts)
    if (size(this%offp_displs) > 0) then
      this%offp_displs(1) = 0
      do i = 2, size(this%offp_displs)
        this%offp_displs(i) = this%offp_displs(i-1) + this%offp_counts(i-1)
      end do
    end if

    !! Communicate the global off-process indices to their owning ranks.
    allocate(this%onp_index(sum(this%onp_counts)))
    call MPI_Neighbor_alltoallv(this%offp_index, this%offp_counts, this%offp_displs, MPI_INTEGER, &
        this%onp_index, this%onp_counts, this%onp_displs, MPI_INTEGER, this%comm, ierr)
    this%onp_index = this%onp_index - this%first + 1  ! map to local indices
    ASSERT(all(this%onp_index >= 1 .and. this%onp_index <= this%onp_size))

  end subroutine add_offp_index

  elemental function global_index(this, n) result(gid)
    class(index_map), intent(in) :: this
    integer, intent(in) :: n
    integer :: gid
    gid = -1
    if (n < 1) return
    if (n <= this%onp_size) then
      gid = this%first + n - 1
    else if (n <= this%local_size) then
      gid = this%offp_index(n-this%onp_size)
    end if
  end function

end module index_map_type
