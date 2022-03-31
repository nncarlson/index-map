    class(index_map), intent(in) :: this
    __DATA_TYPE__, intent(inout) :: onp_data(:)
    __DATA_TYPE__, intent(in), target :: offp_data(:)

    integer :: j, k, n, m

    type box
      __DATA_TYPE__, pointer :: data(:)
    end type
    type(box), allocatable :: offp[:]

    ASSERT(size(onp_data) >= this%onp_size)
    ASSERT(size(offp_data) >= this%offp_size)

    allocate(offp[*])
    offp%data => offp_data
    sync all
    n = 0
    do j = 1, size(this%onp_count)
      associate (i => this%onp_image(j), offset => this%offp_offset(j))
        do k = 1, this%onp_count(j)
          n = n + 1
          m = this%onp_index(n)
          onp_data(m) = __OP__(onp_data(m), offp[i]%data(offset+k))
        end do
      end associate
    end do
#undef __DATA_TYPE__
#undef __OP__
