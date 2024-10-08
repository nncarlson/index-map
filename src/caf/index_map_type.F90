!!
!! INDEX_MAP_TYPE
!!
!! This module provides core capabilities that support the use of distributed
!! arrays in MPI-based SPMD programs through the INDEX_MAP derived type that
!! describes the mapping of an array's index set to processes. The mapping
!! allows for overlap between processes, and provides collective gather and
!! scatter procedures associated with that overlap. The module provides
!! additional procedures for distributing and localizing serial indirect
!! indexing arrays.
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

  use,intrinsic :: iso_fortran_env, only: i4 => int32, i8 => int64, r4 => real32, r8 => real64
  use coarray_collectives, only: co_sum_scan
  implicit none
  private

  type, public :: index_map
    integer :: comm ! PROTECTED
    integer, private :: nproc, rank, root=1
    logical, private :: is_root
    integer :: onp_size=0, offp_size=0, local_size=0  ! PROTECTED
    integer :: global_size=0, first_gid=0, last_gid=0 ! PROTECTED; int64 option?
    ! off-process gather/scatter communication data
    integer, allocatable :: offp_index(:)  ! int64 option?
    integer, allocatable, private :: onp_count(:), onp_image(:), onp_offset(:), onp_index(:)
    integer, allocatable, private :: offp_count(:), offp_image(:), offp_offset(:)
  contains
    generic   :: init => init_dist, init_root, init_dist_offp, init_root_offp, init_ragged
    procedure :: add_offp_index
    procedure :: global_index
    generic :: gather_offp => &
        gath1_i4_1, gath2_i4_1, gath1_i4_2, gath2_i4_2, gath1_i4_3, gath2_i4_3, &
        gath1_r4_1, gath2_r4_1, gath1_r4_2, gath2_r4_2, gath1_r4_3, gath2_r4_3, &
        gath1_r8_1, gath2_r8_1, gath1_r8_2, gath2_r8_2, gath1_r8_3, gath2_r8_3, &
        gath1_c4_1, gath2_c4_1, gath1_c4_2, gath2_c4_2, gath1_c4_3, gath2_c4_3, &
        gath1_c8_1, gath2_c8_1, gath1_c8_2, gath2_c8_2, gath1_c8_3, gath2_c8_3, &
        gath1_dl_1, gath2_dl_1, gath1_dl_2, gath2_dl_2, gath1_dl_3, gath2_dl_3
    generic :: scatter_offp_sum => &
        scat1_sum_i4_1, scat2_sum_i4_1, &
        scat1_sum_r4_1, scat2_sum_r4_1, &
        scat1_sum_r8_1, scat2_sum_r8_1, &
        scat1_sum_c4_1, scat2_sum_c4_1, &
        scat1_sum_c8_1, scat2_sum_c8_1
    generic :: scatter_offp_min => &
        scat1_min_i4_1, scat2_min_i4_1, &
        scat1_min_r4_1, scat2_min_r4_1, &
        scat1_min_r8_1, scat2_min_r8_1
    generic :: scatter_offp_max => &
        scat1_max_i4_1, scat2_max_i4_1, &
        scat1_max_r4_1, scat2_max_r4_1, &
        scat1_max_r8_1, scat2_max_r8_1
    generic :: scatter_offp_or => scat1_or_dl_1, scat2_or_dl_1
    generic :: scatter_offp_and => scat1_and_dl_1, scat2_and_dl_1
    generic :: distribute => &
        dist_i4_1, dist_i8_1, dist_r4_1, dist_r8_1, dist_c4_1, dist_c8_1, dist_dl_1, &
        dist_i4_2, dist_i8_2, dist_r4_2, dist_r8_2, dist_c4_2, dist_c8_2, dist_dl_2, &
        dist_i4_3, dist_i8_3, dist_r4_3, dist_r8_3, dist_c4_3, dist_c8_3, dist_dl_3
    generic :: collate => &
        coll_i4_1, coll_i8_1, coll_r4_1, coll_r8_1, coll_c4_1, coll_c8_1, coll_dl_1, &
        coll_i4_2, coll_i8_2, coll_r4_2, coll_r8_2, coll_c4_2, coll_c8_2, coll_dl_2, &
        coll_i4_3, coll_i8_3, coll_r4_3, coll_r8_3, coll_c4_3, coll_c8_3, coll_dl_3
    generic :: localize_index_array => localize_index_array_serial_1, localize_index_array_serial_2, &
        localize_index_array_dist_1, localize_index_array_dist_2, localize_index_struct_serial
    procedure, private :: init_dist, init_root, init_dist_offp, init_root_offp, init_ragged
    procedure, private :: &
        gath1_i4_1, gath2_i4_1, gath1_i4_2, gath2_i4_2, gath1_i4_3, gath2_i4_3, &
        gath1_r4_1, gath2_r4_1, gath1_r4_2, gath2_r4_2, gath1_r4_3, gath2_r4_3, &
        gath1_r8_1, gath2_r8_1, gath1_r8_2, gath2_r8_2, gath1_r8_3, gath2_r8_3, &
        gath1_c4_1, gath2_c4_1, gath1_c4_2, gath2_c4_2, gath1_c4_3, gath2_c4_3, &
        gath1_c8_1, gath2_c8_1, gath1_c8_2, gath2_c8_2, gath1_c8_3, gath2_c8_3, &
        gath1_dl_1, gath2_dl_1, gath1_dl_2, gath2_dl_2, gath1_dl_3, gath2_dl_3
    procedure, private :: &
        scat1_sum_i4_1, scat2_sum_i4_1, &
        scat1_sum_r4_1, scat2_sum_r4_1, &
        scat1_sum_r8_1, scat2_sum_r8_1, &
        scat1_sum_c4_1, scat2_sum_c4_1, &
        scat1_sum_c8_1, scat2_sum_c8_1
    procedure, private :: &
        scat1_min_i4_1, scat2_min_i4_1, &
        scat1_min_r4_1, scat2_min_r4_1, &
        scat1_min_r8_1, scat2_min_r8_1
    procedure, private :: &
        scat1_max_i4_1, scat2_max_i4_1, &
        scat1_max_r4_1, scat2_max_r4_1, &
        scat1_max_r8_1, scat2_max_r8_1
    procedure, private :: scat1_or_dl_1, scat2_or_dl_1
    procedure, private :: scat1_and_dl_1, scat2_and_dl_1
    procedure, private :: &
        dist_i4_1, dist_i8_1, dist_r4_1, dist_r8_1, dist_c4_1, dist_c8_1, dist_dl_1, &
        dist_i4_2, dist_i8_2, dist_r4_2, dist_r8_2, dist_c4_2, dist_c8_2, dist_dl_2, &
        dist_i4_3, dist_i8_3, dist_r4_3, dist_r8_3, dist_c4_3, dist_c8_3, dist_dl_3
    procedure, private :: &
        coll_i4_1, coll_i8_1, coll_r4_1, coll_r8_1, coll_c4_1, coll_c8_1, coll_dl_1, &
        coll_i4_2, coll_i8_2, coll_r4_2, coll_r8_2, coll_c4_2, coll_c8_2, coll_dl_2, &
        coll_i4_3, coll_i8_3, coll_r4_3, coll_r8_3, coll_c4_3, coll_c8_3, coll_dl_3
    procedure, private :: localize_index_array_serial_1, localize_index_array_serial_2, &
        localize_index_array_dist_1, localize_index_array_dist_2, localize_index_struct_serial
    procedure, private :: add_offp_index_set ! Type bound to workaound gfortran bug
  end type

  interface
    module subroutine gath1_i4_1(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gath1_r4_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gath1_r8_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gath1_c4_1(this, local_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gath1_c8_1(this, local_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine gath1_dl_1(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:)
    end subroutine

    module subroutine gath2_i4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: onp_data(:)
      integer(i4), intent(inout), target :: offp_data(:)
    end subroutine
    module subroutine gath2_r4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: onp_data(:)
      real(r4), intent(inout), target :: offp_data(:)
    end subroutine
    module subroutine gath2_r8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: onp_data(:)
      real(r8), intent(inout), target :: offp_data(:)
    end subroutine
    module subroutine gath2_c4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: onp_data(:)
      complex(r4), intent(inout), target :: offp_data(:)
    end subroutine
    module subroutine gath2_c8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: onp_data(:)
      complex(r8), intent(inout), target :: offp_data(:)
    end subroutine
    module subroutine gath2_dl_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(in) :: onp_data(:)
      logical, intent(inout), target :: offp_data(:)
    end subroutine

    module subroutine gath1_i4_2(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gath1_r4_2(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gath1_r8_2(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gath1_c4_2(this, local_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gath1_c8_2(this, local_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(inout) :: local_data(:,:)
    end subroutine
    module subroutine gath1_dl_2(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:,:)
    end subroutine

    module subroutine gath2_i4_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: onp_data(:,:)
      integer(i4), intent(inout), target :: offp_data(:,:)
    end subroutine
    module subroutine gath2_r4_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: onp_data(:,:)
      real(r4), intent(inout), target :: offp_data(:,:)
    end subroutine
    module subroutine gath2_r8_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: onp_data(:,:)
      real(r8), intent(inout), target :: offp_data(:,:)
    end subroutine
    module subroutine gath2_c4_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: onp_data(:,:)
      complex(r4), intent(inout), target :: offp_data(:,:)
    end subroutine
    module subroutine gath2_c8_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: onp_data(:,:)
      complex(r8), intent(inout), target :: offp_data(:,:)
    end subroutine
    module subroutine gath2_dl_2(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(in) :: onp_data(:,:)
      logical, intent(inout), target :: offp_data(:,:)
    end subroutine

    module subroutine gath1_i4_3(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gath1_r4_3(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gath1_r8_3(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gath1_c4_3(this, local_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gath1_c8_3(this, local_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(inout) :: local_data(:,:,:)
    end subroutine
    module subroutine gath1_dl_3(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:,:,:)
    end subroutine

    module subroutine gath2_i4_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: onp_data(:,:,:)
      integer(i4), intent(inout), target :: offp_data(:,:,:)
    end subroutine
    module subroutine gath2_r4_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: onp_data(:,:,:)
      real(r4), intent(inout), target :: offp_data(:,:,:)
    end subroutine
    module subroutine gath2_r8_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: onp_data(:,:,:)
      real(r8), intent(inout), target :: offp_data(:,:,:)
    end subroutine
    module subroutine gath2_c4_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: onp_data(:,:,:)
      complex(r4), intent(inout), target :: offp_data(:,:,:)
    end subroutine
    module subroutine gath2_c8_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: onp_data(:,:,:)
      complex(r8), intent(inout), target :: offp_data(:,:,:)
    end subroutine
    module subroutine gath2_dl_3(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(in) :: onp_data(:,:,:)
      logical, intent(inout), target :: offp_data(:,:,:)
    end subroutine
  end interface

  interface
    module subroutine scat1_sum_i4_1(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_sum_i4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: onp_data(:)
      integer(i4), intent(in), target :: offp_data(:)
    end subroutine
    module subroutine scat1_sum_r4_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_sum_r4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: onp_data(:)
      real(r4), intent(in), target :: offp_data(:)
    end subroutine
    module subroutine scat1_sum_r8_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_sum_r8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: onp_data(:)
      real(r8), intent(in), target :: offp_data(:)
    end subroutine
    module subroutine scat1_sum_c4_1(this, local_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_sum_c4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r4), intent(inout) :: onp_data(:)
      complex(r4), intent(in), target :: offp_data(:)
    end subroutine
    module subroutine scat1_sum_c8_1(this, local_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_sum_c8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      complex(r8), intent(inout) :: onp_data(:)
      complex(r8), intent(in), target :: offp_data(:)
    end subroutine
  end interface

  interface
    module subroutine scat1_min_i4_1(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_min_i4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: onp_data(:)
      integer(i4), intent(in), target :: offp_data(:)
    end subroutine
    module subroutine scat1_min_r4_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_min_r4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: onp_data(:)
      real(r4), intent(in), target :: offp_data(:)
    end subroutine
    module subroutine scat1_min_r8_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_min_r8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: onp_data(:)
      real(r8), intent(in), target :: offp_data(:)
    end subroutine
  end interface

  interface
    module subroutine scat1_max_i4_1(this, local_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_max_i4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      integer(i4), intent(inout) :: onp_data(:)
      integer(i4), intent(in), target :: offp_data(:)
    end subroutine
    module subroutine scat1_max_r4_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_max_r4_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r4), intent(inout) :: onp_data(:)
      real(r4), intent(in), target :: offp_data(:)
    end subroutine
    module subroutine scat1_max_r8_1(this, local_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_max_r8_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      real(r8), intent(inout) :: onp_data(:)
      real(r8), intent(in), target :: offp_data(:)
    end subroutine
  end interface

  interface
    module subroutine scat1_or_dl_1(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_or_dl_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: onp_data(:)
      logical, intent(in), target :: offp_data(:)
    end subroutine
  end interface

  interface
    module subroutine scat1_and_dl_1(this, local_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: local_data(:)
    end subroutine
    module subroutine scat2_and_dl_1(this, onp_data, offp_data)
      class(index_map), intent(in) :: this
      logical, intent(inout) :: onp_data(:)
      logical, intent(in), target :: offp_data(:)
    end subroutine
  end interface

  interface
    module subroutine dist_i4_1(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in), target :: src(:)
      integer(i4), intent(inout) :: dest(:)
    end subroutine
    module subroutine dist_i8_1(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in), target :: src(:)
      integer(i8), intent(inout) :: dest(:)
    end subroutine
    module subroutine dist_r4_1(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in), target :: src(:)
      real(r4), intent(inout) :: dest(:)
    end subroutine
    module subroutine dist_r8_1(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in), target :: src(:)
      real(r8), intent(inout) :: dest(:)
    end subroutine
    module subroutine dist_c4_1(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in), target :: src(:)
      complex(r4), intent(inout) :: dest(:)
    end subroutine
    module subroutine dist_c8_1(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in), target :: src(:)
      complex(r8), intent(inout) :: dest(:)
    end subroutine
    module subroutine dist_dl_1(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in), target :: src(:)
      logical, intent(inout) :: dest(:)
    end subroutine
    module subroutine dist_i4_2(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in), target :: src(:,:)
      integer(i4), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine dist_i8_2(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in), target :: src(:,:)
      integer(i8), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine dist_r4_2(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in), target :: src(:,:)
      real(r4), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine dist_r8_2(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in), target :: src(:,:)
      real(r8), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine dist_c4_2(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in), target :: src(:,:)
      complex(r4), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine dist_c8_2(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in), target :: src(:,:)
      complex(r8), intent(inout) :: dest(:,:)
    end subroutine
    module subroutine dist_dl_2(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in), target :: src(:,:)
      logical, intent(inout) :: dest(:,:)
    end subroutine
    module subroutine dist_i4_3(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in), target :: src(:,:,:)
      integer(i4), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine dist_i8_3(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in), target :: src(:,:,:)
      integer(i8), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine dist_r4_3(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in), target :: src(:,:,:)
      real(r4), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine dist_r8_3(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in), target :: src(:,:,:)
      real(r8), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine dist_c4_3(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in), target :: src(:,:,:)
      complex(r4), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine dist_c8_3(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in), target :: src(:,:,:)
      complex(r8), intent(inout) :: dest(:,:,:)
    end subroutine
    module subroutine dist_dl_3(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in), target :: src(:,:,:)
      logical, intent(inout) :: dest(:,:,:)
    end subroutine
  end interface

  interface
    module subroutine coll_i4_1(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: src(:)
      integer(i4), intent(inout), target :: dest(:)
    end subroutine
    module subroutine coll_i8_1(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in) :: src(:)
      integer(i8), intent(inout), target :: dest(:)
    end subroutine
    module subroutine coll_r4_1(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: src(:)
      real(r4), intent(inout), target :: dest(:)
    end subroutine
    module subroutine coll_r8_1(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: src(:)
      real(r8), intent(inout), target :: dest(:)
    end subroutine
    module subroutine coll_c4_1(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: src(:)
      complex(r4), intent(inout), target :: dest(:)
    end subroutine
    module subroutine coll_c8_1(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: src(:)
      complex(r8), intent(inout),target :: dest(:)
    end subroutine
    module subroutine coll_dl_1(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in) :: src(:)
      logical, intent(inout), target :: dest(:)
    end subroutine
    module subroutine coll_i4_2(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: src(:,:)
      integer(i4), intent(inout), target :: dest(:,:)
    end subroutine
    module subroutine coll_i8_2(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in) :: src(:,:)
      integer(i8), intent(inout), target :: dest(:,:)
    end subroutine
    module subroutine coll_r4_2(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: src(:,:)
      real(r4), intent(inout), target :: dest(:,:)
    end subroutine
    module subroutine coll_r8_2(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: src(:,:)
      real(r8), intent(inout), target :: dest(:,:)
    end subroutine
    module subroutine coll_c4_2(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: src(:,:)
      complex(r4), intent(inout), target :: dest(:,:)
    end subroutine
    module subroutine coll_c8_2(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: src(:,:)
      complex(r8), intent(inout), target :: dest(:,:)
    end subroutine
    module subroutine coll_dl_2(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in) :: src(:,:)
      logical, intent(inout), target :: dest(:,:)
    end subroutine
    module subroutine coll_i4_3(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i4), intent(in) :: src(:,:,:)
      integer(i4), intent(inout), target :: dest(:,:,:)
    end subroutine
    module subroutine coll_i8_3(this, src, dest)
      class(index_map), intent(in) :: this
      integer(i8), intent(in) :: src(:,:,:)
      integer(i8), intent(inout), target :: dest(:,:,:)
    end subroutine
    module subroutine coll_r4_3(this, src, dest)
      class(index_map), intent(in) :: this
      real(r4), intent(in) :: src(:,:,:)
      real(r4), intent(inout), target :: dest(:,:,:)
    end subroutine
    module subroutine coll_r8_3(this, src, dest)
      class(index_map), intent(in) :: this
      real(r8), intent(in) :: src(:,:,:)
      real(r8), intent(inout), target :: dest(:,:,:)
    end subroutine
    module subroutine coll_c4_3(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r4), intent(in) :: src(:,:,:)
      complex(r4), intent(inout), target :: dest(:,:,:)
    end subroutine
    module subroutine coll_c8_3(this, src, dest)
      class(index_map), intent(in) :: this
      complex(r8), intent(in) :: src(:,:,:)
      complex(r8), intent(inout), target :: dest(:,:,:)
    end subroutine
    module subroutine coll_dl_3(this, src, dest)
      class(index_map), intent(in) :: this
      logical, intent(in) :: src(:,:,:)
      logical, intent(inout), target :: dest(:,:,:)
    end subroutine
  end interface

  interface
    module subroutine localize_index_array_serial_1(domain, g_index, range, l_index, stat)
      class(index_map), intent(in) :: domain
      class(index_map), intent(inout) :: range
      integer, intent(in) :: g_index(:)
      integer, allocatable, intent(out) :: l_index(:)
      integer, intent(out), optional :: stat
    end subroutine
    module subroutine localize_index_array_serial_2(domain, g_index, range, l_index, stat)
      class(index_map), intent(in) :: domain
      class(index_map), intent(inout) :: range
      integer, intent(in) :: g_index(:,:)
      integer, allocatable, intent(out) :: l_index(:,:)
      integer, intent(out), optional :: stat
    end subroutine
    module subroutine localize_index_array_dist_1(range, index, stat)
      class(index_map), intent(inout) :: range
      integer, intent(inout) :: index(:)
      integer, intent(out), optional :: stat
    end subroutine
    module subroutine localize_index_array_dist_2(range, index, stat)
      class(index_map), intent(inout) :: range
      integer, contiguous, intent(inout), target :: index(:,:)
      integer, intent(out), optional :: stat
    end subroutine
    module subroutine localize_index_struct_serial(domain, g_count, g_index, range, l_count, l_index, stat)
      class(index_map), intent(in) :: domain
      class(index_map), intent(inout) :: range
      integer, intent(in) :: g_index(:), g_count(:)
      integer, allocatable, intent(out) :: l_index(:), l_count(:)
      integer, intent(out), optional :: stat
    end subroutine
  end interface

contains

  !! Each image supplies its block size
  subroutine init_dist(this, bsize, root)

    class(index_map), intent(out) :: this
    integer, intent(in) :: bsize
    integer, intent(in), optional :: root

    this%nproc = num_images()
    this%rank = this_image()
    if (present(root)) this%root = root ! overwrite default
    this%is_root = (this%rank == this%root)

    this%onp_size = bsize
    this%offp_size = 0
    this%local_size = this%onp_size + this%offp_size
    call co_sum_scan(bsize, this%last_gid)
    this%first_gid = this%last_gid - this%onp_size + 1
    this%global_size = this%last_gid
    call co_broadcast(this%global_size, this%nproc)

  end subroutine

  !! One root rank has an array of image block sizes
  subroutine init_root(this, bsizes, root)
    class(index_map), intent(out) :: this
    integer, intent(in) :: bsizes(:)
    integer, intent(in), optional :: root
    integer :: i
    integer, allocatable :: bsize[:]
    allocate(bsize[*])
    if (present(root)) this%root = root ! overwrite default
    if (this_image() == this%root) then
      INSIST(size(bsizes) == num_images())
      do i = 1, num_images()
        bsize[i] = bsizes(i)
      end do
    end if
    sync all
    call init_dist(this, bsize, root)
  end subroutine

  !! Each rank supplied with its block size and list of off-process indices
  subroutine init_dist_offp(this, bsize, offp_index, root)
    class(index_map), intent(out) :: this
    integer, intent(in) :: bsize
    integer, intent(in) :: offp_index(:) ! 64-bit option?
    integer, intent(in), optional :: root
    call init_dist(this, bsize, root)
    call add_offp_index(this, offp_index)
  end subroutine

  !! One root rank has an array of rank block sizes and off-process indices for all ranks
  subroutine init_root_offp(this, bsizes, offp_counts, offp_indices, root)
    class(index_map), intent(out) :: this
    integer, intent(in) :: bsizes(:), offp_counts(:)
    integer, intent(in) :: offp_indices(:)
    integer, intent(in), optional :: root
    type(index_map) :: temp
    integer, allocatable :: offp_index(:)
    call init_root(this, bsizes, root)
    call init_root(temp, offp_counts, root)
    allocate(offp_index(temp%onp_size))
    call temp%distribute(offp_indices, offp_index)
    call add_offp_index(this, offp_index)
  end subroutine

  !! Ragged index set based on another index set
  subroutine init_ragged(this, domain, g_count)

    class(index_map), intent(out) :: this
    type(index_map), intent(in) :: domain
    integer, intent(in) :: g_count(:)

    integer :: n, i, j, nmax
    integer, allocatable :: l_count(:), offset(:), offp_index(:)

    ASSERT(size(g_count) >= merge(domain%global_size,0,domain%is_root))

    allocate(l_count(domain%local_size))
    call domain%distribute(g_count, l_count)
    if (allocated(domain%offp_index)) call domain%gather_offp(l_count)

    call this%init(sum(l_count(1:domain%onp_size)))

    if (allocated(domain%offp_index)) then
      call co_sum_scan(this%onp_size, n)
      n = n - this%onp_size ! exclusive scan of this%onp_size
      allocate(offset(domain%local_size))
      if (size(offset) > 0) then
        offset(1) = n
        do j = 2, domain%onp_size
          offset(j) = offset(j-1) + l_count(j-1)
        end do
        call domain%gather_offp(offset)
      end if
      n = sum(l_count(domain%onp_size+1:domain%local_size))
      nmax = n; call co_max(nmax)
      if (nmax > 0) then ! have off-process indices for this index map
        allocate(offp_index(n))
        n = 0
        do j = domain%onp_size + 1, domain%local_size
          offp_index(n+1:n+l_count(j)) = [(offset(j)+i, i=1,l_count(j))]
          n = n + l_count(j)
        end do
        call add_offp_index(this, offp_index)
      end if
    end if

  end subroutine

  !! Return the global index corresponding to local index N
  elemental function global_index(this, n) result(gid)
    class(index_map), intent(in) :: this
    integer, intent(in) :: n
    integer :: gid
    gid = -1
    if (n < 1) return
    if (n <= this%onp_size) then
      gid = this%first_gid + n - 1
    else if (n <= this%local_size) then
      gid = this%offp_index(n-this%onp_size)
    end if
  end function

  !! Add off-process indices to an already initialized index map.
  !! This is a public interface with untrusted OFFP_INDEX input.

  subroutine add_offp_index(this, offp_index)
    use integer_set_type
    class(index_map), intent(inout) :: this
    integer, intent(in) :: offp_index(:)
    type(integer_set) :: offp_set
    ASSERT(minval(offp_index) >= 1)
    ASSERT(maxval(offp_index) <= this%global_size)
    ASSERT(all((offp_index < this%first_gid) .or. (offp_index > this%last_gid)))
    call offp_set%add(offp_index) ! sort and remove duplicates
    call add_offp_index_set(this, offp_set)
  end subroutine

  !! Add off-process indices to an already initialized index map.
  !! This is an internal interface with trusted OFFP_SET input.
  !!
  !! This defines the following components which are used by the GATHER_OFFP
  !! and SCATTER_OFFP type bound subroutines.
  !!
  !! OFFP_INDEX is the array of off-process global indices included on this
  !! image along with the on-process indices assigned to the image. The array
  !! is ordered and has an associated block partitioning described by the pair
  !! of conformable arrays OFFP_IMAGE and OFFP_COUNT: the initial block of
  !! OFFP_COUNT(1) elements are owned by image OFFP_IMAGE(1), the next block
  !! of OFFP_COUNT(2) elements are owned by OFFP_IMAGE(2), and so forth. The
  !! values of OFFP_IMAGE are strictly increasing.
  !!
  !! ONP_INDEX is the array of local on-process indices that are present on
  !! other images as off-process indices. The array has an associated block
  !! partitioning described by the pair of conformable arrays ONP_IMAGE and
  !! ONP_COUNT: the initial block of ONP_COUNT(1) elements are included as
  !! off-process elements on image ONP_IMAGE(1), the next block of
  !! ONP_COUNT(2) elements are included as off-process on image ONP_IMAGE(2),
  !! and so forth. The values of ONP_IMAGE are strictly increasing.
  !!
  !! The blocks of the ONP_INDEX arrays correspond, one-to-one, to blocks of
  !! the OFFP_INDEX arrays. On image i let i' = ONP_IMAGE(j) for some block j.
  !! Then on image i' there is a unique j' such that i = OFFP_IMAGE(j'). And
  !! conversely, if on image i', i = ONP_IMAGE(j') for some block j', then
  !! on image i there is a unique j such that i' = OFFP_IMAGE(j). And then
  !! block j of ONP_INDEX on image i corresponds, element for element, to
  !! block j' of OFFP_INDEX on image i'.
  !!
  !! The arrays OFFP_OFFSET and ONP_OFFSET facilitate communication between
  !! corresponding blocks. OFFP_OFFSET is conformable with ONP_COUNT and
  !! ONP_IMAGE. On image i, OFFP_OFFSET(j) is the offset into OFFP_INDEX
  !! on image i' = ONP_IMAGE(j) for the start of the remote block that
  !! corresponds to ONP_INDEX block j. Similiarly, ONP_OFFSET is conformable
  !! with OFFP_COUNT and OFFP_IMAGE. On image i', ONP_OFFSET(j') is the
  !! offset into ONP_INDEX on image i = OFFP_IMAGE(j') for the start of the
  !! remote block that corresponds to OFFP_INDEX block j'.

  subroutine add_offp_index_set(this, offp_set)

    use integer_set_type

    class(index_map), intent(inout), target :: this
    type(integer_set), intent(inout) :: offp_set

    integer :: i, j, j1, k, n
    integer, allocatable :: offp_image(:), onp_counts(:)
    integer, allocatable :: last[:], offp_counts(:)[:], offsets(:)[:]

    type box
      integer, pointer :: data(:)
    end type
    type(box), allocatable :: buffer[:]

    !TODO? Allow extending an existing %offp_index
    INSIST(.not.allocated(this%offp_index))

    this%offp_index = offp_set
    call offp_set%clear
    this%offp_size  = size(this%offp_index)
    this%local_size = this%onp_size + this%offp_size

    !! Determine the image that owns each off-process index (OFFP_IMAGE).
    !! NB: This assumes OFFP_INDEX is ordered; OFFP_IMAGE will be ordered.
    allocate(last[*], offp_image(this%offp_size))
    last = this%last_gid
    sync all
    i = 1
    do j = 1, size(this%offp_index)
      do while (this%offp_index(j) > last[i])
        i = i + 1
        INSIST(i <= this%nproc)
      end do
      INSIST(i /= this_image())
      offp_image(j) = i
    end do
    
    !! Get the number of off-process indices owned by each image (OFFP_COUNTS).
    !! NB: This assumes the OFFP_IMAGE array is ordered.
    allocate(offp_counts(this%nproc)[*])
    j1 = 1
    do i = 1, this%nproc
      do j = j1, size(offp_image)
        if (offp_image(j) > i) exit
      end do
      offp_counts(i) = j - j1
      j1 = j
    end do

    !! Communicate the count of off-process indices owned by an image to that
    !! image. The result ONP_COUNTS is the number of on-process indices that
    !! are included as off-process indices on each of the other images.
    sync all
    allocate(onp_counts(this%nproc))
    do i = 1, this%nproc
      onp_counts(i) = offp_counts(this_image())[i]
    end do
    
    !! Offsets into the ONP_INDEX array for each of the image blocks.
    allocate(offsets(this%nproc)[*])
    offsets(1) = 0
    do i = 2, this%nproc
      offsets(i) = offsets(i-1) + onp_counts(i-1)
    end do
    
    !! Compress the OFFP_COUNTS array, dropping elements with 0 count.
    sync all
    n = count(offp_counts > 0)
    allocate(this%offp_count(n), this%offp_image(n), this%onp_offset(n))
    n = 0
    do i = 1, this%nproc
      if (offp_counts(i) == 0) cycle
      n = n + 1
      this%offp_image(n) = i
      this%offp_count(n) = offp_counts(i)
      this%onp_offset(n) = offsets(this_image())[i]
    end do
    
    !! Offsets into the OFFP_INDEX array for each of the image blocks.
    sync all
    offsets(1) = 0
    do i = 2, this%nproc
      offsets(i) = offsets(i-1) + offp_counts(i-1)
    end do

    !! Compress the ONP_COUNTS array, dropping elements with 0 count.
    sync all
    n = count(onp_counts > 0)
    allocate(this%onp_count(n), this%onp_image(n), this%offp_offset(n))
    n = 0
    do i = 1, this%nproc
      if (onp_counts(i) == 0) cycle
      n = n + 1
      this%onp_image(n) = i
      this%onp_count(n) = onp_counts(i)
      this%offp_offset(n) = offsets(this_image())[i]
    end do
    
    !! Communicate the global off-process indices to their owning images
    !! and map to the corresponding local on-process indices (ONP_INDEX).
    n = sum(this%onp_count)
    allocate(this%onp_index(n), buffer[*])
    buffer%data => this%offp_index
    sync all
    n = 0
    do j = 1, size(this%onp_count)
      associate (i => this%onp_image(j), offset => this%offp_offset(j))
        do k = 1, this%onp_count(j)
          n = n + 1
          this%onp_index(n) = buffer[i]%data(offset+k) - this%first_gid + 1
        end do
      end associate
    end do
    ASSERT(all(this%onp_index >= 1 .and. this%onp_index <= this%onp_size))

    ASSERT(gather_offp_verified(this))

  end subroutine add_offp_index_set

  !! This function returns true if the gather_offp operation returns the
  !! expected result for integer data, and otherwise it returns false.
  !! Useful for testing with live index_map data.

  logical function gather_offp_verified(this) result(pass)
    type(index_map), intent(in) :: this
    integer :: j, n, onp_data(this%onp_size), offp_data(this%offp_size)
    do j = 1, this%onp_size
      onp_data(j) = global_index(this, j)
    end do
    call this%gather_offp(onp_data, offp_data)
    n = merge(1, 0, all(offp_data == this%offp_index))
    call co_min(n)
    pass = (n == 1)
  end function

end module index_map_type
