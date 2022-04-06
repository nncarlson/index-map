    use,intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
    class(index_map), intent(inout) :: this
    __DATA_TYPE__, intent(in) :: onp_data(:)
    __DATA_TYPE__, intent(inout), target :: offp_data(:)

    integer :: j, k, n
    __DATA_TYPE__, pointer :: onp_buffer(:)

    ASSERT(size(onp_data) >= this%onp_size)
    ASSERT(size(offp_data) >= this%offp_size)

    if (.not.allocated(this%offp_index)) return

    call c_f_pointer(c_loc(this%onp_buffer(1)), onp_buffer, shape=[size(this%onp_index)])
    
    do j = 1, size(this%onp_index)
      onp_buffer(j) = onp_data(this%onp_index(j))
    end do

    co_view%data => offp_data
    sync all

    n = 0
    do j = 1, size(this%onp_count)
      associate (i => this%onp_image(j), offset => this%offp_offset(j), m => this%onp_count(j))
        co_view[i]%data(offset+1:offset+m) = onp_buffer(n+1:n+m)
        n = n + m
      end associate
    end do
    sync all
#undef __DATA_TYPE__
