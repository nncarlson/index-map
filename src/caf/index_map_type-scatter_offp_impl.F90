!! Implementation of INDEX_MAP SCATTER_OFFP Procedures
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

submodule(index_map_type) scatter_impl
implicit none
contains

  module subroutine scat1_sum_i4_1(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:)
    call scat2_sum_i4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_sum_i4_1(this, onp_data, offp_data)
#define __DATA_TYPE__ integer(i4)
#define __OP__(a,b) a + b
#include "scat_1.inc"
  end subroutine

  module subroutine scat1_sum_r4_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:)
    call scat2_sum_r4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_sum_r4_1(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r4)
#define __OP__(a,b) a + b
#include "scat_1.inc"
  end subroutine

  module subroutine scat1_sum_r8_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:)
    call scat2_sum_r8_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_sum_r8_1(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r8)
#define __OP__(a,b) a + b
#include "scat_1.inc"
  end subroutine

  module subroutine scat1_min_i4_1(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:)
    call scat2_min_i4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_min_i4_1(this, onp_data, offp_data)
#define __DATA_TYPE__ integer(i4)
#define __OP__(a,b) min(a, b)
#include "scat_1.inc"
  end subroutine

  module subroutine scat1_min_r4_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:)
    call scat2_min_r4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_min_r4_1(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r4)
#define __OP__(a,b) min(a, b)
#include "scat_1.inc"
  end subroutine

  module subroutine scat1_min_r8_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:)
    call scat2_min_r8_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_min_r8_1(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r8)
#define __OP__(a,b) min(a, b)
#include "scat_1.inc"
  end subroutine

  module subroutine scat1_max_i4_1(this, local_data)
    class(index_map), intent(in) :: this
    integer(i4), intent(inout) :: local_data(:)
    call scat2_max_i4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_max_i4_1(this, onp_data, offp_data)
#define __DATA_TYPE__ integer(i4)
#define __OP__(a,b) max(a, b)
#include "scat_1.inc"
  end subroutine

  module subroutine scat1_max_r4_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r4), intent(inout) :: local_data(:)
    call scat2_max_r4_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_max_r4_1(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r4)
#define __OP__(a,b) max(a, b)
#include "scat_1.inc"
  end subroutine

  module subroutine scat1_max_r8_1(this, local_data)
    class(index_map), intent(in) :: this
    real(r8), intent(inout) :: local_data(:)
    call scat2_max_r8_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_max_r8_1(this, onp_data, offp_data)
#define __DATA_TYPE__ real(r8)
#define __OP__(a,b) max(a, b)
#include "scat_1.inc"
  end subroutine

  module subroutine scat1_or_dl_1(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:)
    call scat2_or_dl_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_or_dl_1(this, onp_data, offp_data)
#define __DATA_TYPE__ logical
#define __OP__(a,b) a .or. b
#include "scat_1.inc"
  end subroutine

  module subroutine scat1_and_dl_1(this, local_data)
    class(index_map), intent(in) :: this
    logical, intent(inout) :: local_data(:)
    call scat2_and_dl_1(this, local_data(:this%onp_size), local_data(this%onp_size+1:))
  end subroutine

  module subroutine scat2_and_dl_1(this, onp_data, offp_data)
#define __DATA_TYPE__ logical
#define __OP__(a,b) a .and. b
#include "scat_1.inc"
  end subroutine

end submodule scatter_impl
