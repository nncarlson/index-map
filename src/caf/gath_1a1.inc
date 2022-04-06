    class(index_map), intent(inout) :: this
    __DATA_TYPE__, intent(in) :: onp_data(:)
    __DATA_TYPE__, intent(inout), target :: offp_data(:)

    integer :: j, k, n

    ASSERT(size(onp_data) >= this%onp_size)
    ASSERT(size(offp_data) >= this%offp_size)

    if (.not.allocated(this%offp_index)) return

    this%co_buffer%data => offp_data
    sync all

    n = 0
    do j = 1, size(this%onp_count)
      associate (i => this%onp_image(j), offset => this%offp_offset(j))
        do k = 1, this%onp_count(j)
          n = n + 1
          this%co_buffer[i]%data(offset+k) = onp_data(this%onp_index(n))
        end do
      end associate
    end do
    sync all
#undef __DATA_TYPE__
