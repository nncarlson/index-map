    class(index_map), intent(in) :: this
    __DATA_TYPE__, intent(in), target :: onp_data(:)
    __DATA_TYPE__, intent(inout) :: offp_data(:)

    integer :: j, k, n

    ASSERT(size(onp_data) >= this%onp_size)
    ASSERT(size(offp_data) >= this%offp_size)

    if (.not.allocated(this%offp_index)) return

    co_view%data => onp_data
    sync all

    n = 0
    do j = 1, size(this%offp_count)
      associate (i => this%offp_image(j))
        do k = 1, this%offp_count(j)
          n = n + 1
          offp_data(n) = co_view[i]%data(this%src_index(n))
        end do
      end associate
    end do
#undef __DATA_TYPE__
