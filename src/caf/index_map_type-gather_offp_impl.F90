!! Implementation of INDEX_MAP GATHER_OFFP Procedures
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

submodule(index_map_type) gather_offp_impl
implicit none
contains

!!!! I4 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath1_i4_1(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:)
    call gath2_i4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gath2_i4_1(this, onp_data, offp_data)
#define __DATA_TYPE__ integer(i4)
#ifdef __GFORTRAN__
#include "gath_1_alt.inc"
#else
#include "gath_1.inc"
#endif
  end subroutine

  module subroutine gath1_i4_2(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:,:)
    call gath2_i4_2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_i4_2(this, onp_data, offp_data)
#define __DATA_TYPE__ integer(i4)
#include "gath_2.inc"
  end subroutine

  module subroutine gath1_i4_3(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:,:,:)
    call gath2_i4_3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_i4_3(this, onp_data, offp_data)
#define __DATA_TYPE__ integer(i4)
#include "gath_3.inc"
  end subroutine

!!!! R4 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath1_r4_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:)
    call gath2_r4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gath2_r4_1(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r4)
#ifdef __GFORTRAN__
#include "gath_1_alt.inc"
#else
#include "gath_1.inc"
#endif
  end subroutine

  module subroutine gath1_r4_2(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:,:)
    call gath2_r4_2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_r4_2(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r4)
#include "gath_2.inc"
  end subroutine

  module subroutine gath1_r4_3(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:,:,:)
    call gath2_r4_3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_r4_3(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r4)
#include "gath_3.inc"
  end subroutine

!!!! R8 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath1_r8_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:)
    call gath2_r8_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gath2_r8_1(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r8)
#ifdef __GFORTRAN__
#include "gath_1_alt.inc"
#else
#include "gath_1.inc"
#endif
  end subroutine

  module subroutine gath1_r8_2(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:,:)
    call gath2_r8_2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_r8_2(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r8)
#include "gath_2.inc"
  end subroutine

  module subroutine gath1_r8_3(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:,:,:)
    call gath2_r8_3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_r8_3(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r8)
#include "gath_3.inc"
  end subroutine

!!!! C4 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath1_c4_1(this, local_data)
    class(index_map), intent(in) :: this
    complex(r4), intent(inout) :: local_data(:)
    call gath2_c4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gath2_c4_1(this, onp_data, offp_data)
#define __DATA_TYPE__ complex(r4)
#ifdef __GFORTRAN__
#include "gath_1_alt.inc"
#else
#include "gath_1.inc"
#endif
  end subroutine

  module subroutine gath1_c4_2(this, local_data)
    class(index_map), intent(in) :: this
    complex(r4), intent(inout) :: local_data(:,:)
    call gath2_c4_2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_c4_2(this, onp_data, offp_data)
#define __DATA_TYPE__ complex(r4)
#include "gath_2.inc"
  end subroutine

  module subroutine gath1_c4_3(this, local_data)
    class(index_map), intent(in) :: this
    complex(r4), intent(inout) :: local_data(:,:,:)
    call gath2_c4_3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_c4_3(this, onp_data, offp_data)
#define __DATA_TYPE__ complex(r4)
#include "gath_3.inc"
  end subroutine

!!!! C8 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath1_c8_1(this, local_data)
    class(index_map), intent(in) :: this
    complex(r8), intent(inout) :: local_data(:)
    call gath2_c8_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gath2_c8_1(this, onp_data, offp_data)
#define __DATA_TYPE__ complex(r8)
#ifdef __GFORTRAN__
#include "gath_1_alt.inc"
#else
#include "gath_1.inc"
#endif
  end subroutine

  module subroutine gath1_c8_2(this, local_data)
    class(index_map), intent(in) :: this
    complex(r8), intent(inout) :: local_data(:,:)
    call gath2_c8_2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_c8_2(this, onp_data, offp_data)
#define __DATA_TYPE__ complex(r8)
#include "gath_2.inc"
  end subroutine

  module subroutine gath1_c8_3(this, local_data)
    class(index_map), intent(in) :: this
    complex(r8), intent(inout) :: local_data(:,:,:)
    call gath2_c8_3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_c8_3(this, onp_data, offp_data)
#define __DATA_TYPE__ complex(r8)
#include "gath_3.inc"
  end subroutine

!!!! LOGICAL DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine gath1_dl_1(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:)
    call gath2_dl_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine gath2_dl_1(this, onp_data, offp_data)
#define __DATA_TYPE__ logical
#ifdef __GFORTRAN__
#include "gath_1_alt.inc"
#else
#include "gath_1.inc"
#endif
  end subroutine

  module subroutine gath1_dl_2(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:,:)
    call gath2_dl_2(this, local_data(:,:this%onp_size), local_data(:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_dl_2(this, onp_data, offp_data)
#define __DATA_TYPE__ logical
#include "gath_2.inc"
  end subroutine

  module subroutine gath1_dl_3(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:,:,:)
    call gath2_dl_3(this, local_data(:,:,:this%onp_size), local_data(:,:,this%onp_size+1:))
  end subroutine

  module subroutine gath2_dl_3(this, onp_data, offp_data)
#define __DATA_TYPE__ logical
#include "gath_3.inc"
  end subroutine

end submodule gather_offp_impl
