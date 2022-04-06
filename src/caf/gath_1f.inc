    class(index_map), intent(in) :: this
    __DATA_TYPE__, intent(in), target :: onp_data(:)
    __DATA_TYPE__, intent(inout) :: offp_data(:)

    integer :: j, k, n

    ASSERT(size(onp_data) >= this%onp_size)
    ASSERT(size(offp_data) >= this%offp_size)

    if (.not.allocated(this%offp_index)) return

    do j = 1, size(this%onp_index)
      co_onp%buffer(j) = onp_data(this%onp_index(j))
    end do
    sync all

    n = 0
    do j = 1, size(this%offp_count)
      associate (i => this%offp_image(j), offset => this%onp_offset(j), m => this%offp_count(j))
        offp_data(n+1:n+m) = co_onp[i]%buffer(offset+1:offset+m)
        n = n + m
      end associate
    end do
#undef __DATA_TYPE__
