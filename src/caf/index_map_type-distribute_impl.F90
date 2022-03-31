!! Implementation of INDEX_MAP DISTRIBUTE Subroutines
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

submodule(index_map_type) distribute_impl
implicit none
contains

!!!! RANK-1 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine dist_i4_1(this, src, dest)
#define __DATA_TYPE__ integer(i4)
#include "dist_1.inc"
  end subroutine

  module subroutine dist_i8_1(this, src, dest)
#define __DATA_TYPE__ integer(i8)
#include "dist_1.inc"
  end subroutine

  module subroutine dist_r4_1(this, src, dest)
#define __DATA_TYPE__ real(r4)
#include "dist_1.inc"
  end subroutine

  module subroutine dist_r8_1(this, src, dest)
#define __DATA_TYPE__ real(r8)
#include "dist_1.inc"
  end subroutine

  module subroutine dist_dl_1(this, src, dest)
#define __DATA_TYPE__ logical
#include "dist_1.inc"
  end subroutine

!!!! RANK-2 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine dist_i4_2(this, src, dest)
#define __DATA_TYPE__ integer(i4)
#include "dist_2.inc"
  end subroutine

  module subroutine dist_i8_2(this, src, dest)
#define __DATA_TYPE__ integer(i8)
#include "dist_2.inc"
  end subroutine

  module subroutine dist_r4_2(this, src, dest)
#define __DATA_TYPE__ real(r4)
#include "dist_2.inc"
  end subroutine

  module subroutine dist_r8_2(this, src, dest)
#define __DATA_TYPE__ real(r8)
#include "dist_2.inc"
  end subroutine

  module subroutine dist_dl_2(this, src, dest)
#define __DATA_TYPE__ logical
#include "dist_2.inc"
  end subroutine

!!!! RANK-3 DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  module subroutine dist_i4_3(this, src, dest)
#define __DATA_TYPE__ integer(i4)
#include "dist_3.inc"
  end subroutine

  module subroutine dist_i8_3(this, src, dest)
#define __DATA_TYPE__ integer(i8)
#include "dist_3.inc"
  end subroutine

  module subroutine dist_r4_3(this, src, dest)
#define __DATA_TYPE__ real(r4)
#include "dist_3.inc"
  end subroutine

  module subroutine dist_r8_3(this, src, dest)
#define __DATA_TYPE__ real(r8)
#include "dist_3.inc"
  end subroutine

  module subroutine dist_dl_3(this, src, dest)
#define __DATA_TYPE__ logical
#include "dist_3.inc"
  end subroutine

end submodule
