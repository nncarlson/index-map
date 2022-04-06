    class(index_map), intent(in) :: this
    __DATA_TYPE__, intent(in), target :: onp_data(:)
    __DATA_TYPE__, intent(inout) :: offp_data(:)

    integer :: j, k, n

    type box
      __DATA_TYPE__, allocatable :: data(:)
    end type
    type(box), allocatable :: onp[:]

    ASSERT(size(onp_data) >= this%onp_size)
    ASSERT(size(offp_data) >= this%offp_size)

    if (.not.allocated(this%offp_index)) return

    allocate(onp[*])
    allocate(onp%data(size(this%onp_index)))
    do j = 1, size(this%onp_index)
      onp%data(j) = onp_data(this%onp_index(j))
    end do
    sync all

    n = 0
    do j = 1, size(this%offp_count)
      associate (i => this%offp_image(j), offset => this%onp_offset(j), m => this%offp_count(j))
        offp_data(n+1:n+m) = onp[i]%data(offset+1:offset+m)
        n = n + m
      end associate
    end do
#ifdef __GFORTRAN__
    deallocate(onp)
#endif
#undef __DATA_TYPE__
