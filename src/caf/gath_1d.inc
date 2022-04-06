    class(index_map), intent(in) :: this
    __DATA_TYPE__, intent(in), target :: onp_data(:)
    __DATA_TYPE__, intent(inout) :: offp_data(:)

    integer :: j

    ASSERT(size(onp_data) >= this%onp_size)
    ASSERT(size(offp_data) >= this%offp_size)

    if (.not.allocated(this%offp_index)) return

    co_view%data => onp_data
    sync all

    do j = 1, this%offp_size
      offp_data(j) = co_view[this%src_image(j)]%data(this%src_index(j))
    end do
#undef __DATA_TYPE__
