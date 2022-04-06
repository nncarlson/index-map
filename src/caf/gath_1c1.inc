    class(index_map), intent(in) :: this
    __DATA_TYPE__, intent(in) :: onp_data(:)
    __DATA_TYPE__, intent(inout), target :: offp_data(:)

    integer :: j, k, n
    __DATA_TYPE__, allocatable :: onp_buff(:)

    ASSERT(size(onp_data) >= this%onp_size)
    ASSERT(size(offp_data) >= this%offp_size)

    if (.not.allocated(this%offp_index)) return

    allocate(onp_buff(size(this%onp_index)))
    do j = 1, size(this%onp_index)
      onp_buff(j) = onp_data(this%onp_index(j))
    end do

    co_view%data => offp_data
    sync all

    n = 0
    do j = 1, size(this%onp_count)
      associate (i => this%onp_image(j), offset => this%offp_offset(j), m => this%onp_count(j))
        co_view[i]%data(offset+1:offset+m) = onp_buff(n+1:n+m)
        n = n + m
      end associate
    end do
    sync all
#undef __DATA_TYPE__
