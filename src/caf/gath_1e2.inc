    class(index_map), intent(in) :: this
    __DATA_TYPE__, intent(in), target :: onp_data(:)
    __DATA_TYPE__, intent(inout) :: offp_data(:)

    integer :: j, k, n

    type box
      __DATA_TYPE__, pointer :: data(:)
    end type
    type(box), allocatable :: onp[:]
    __DATA_TYPE__, target :: buffer(size(this%onp_index))

    ASSERT(size(onp_data) >= this%onp_size)
    ASSERT(size(offp_data) >= this%offp_size)

    if (.not.allocated(this%offp_index)) return

    !allocate(buffer(size(this%onp_index)))
    do j = 1, size(this%onp_index)
      buffer(j) = onp_data(this%onp_index(j))
    end do
    allocate(onp[*])
    onp%data => buffer
    sync all

    n = 0
    do j = 1, size(this%offp_count)
      associate (i => this%offp_image(j), offset => this%onp_offset(j))
        do k = 1, this%offp_count(j)
          n = n + 1
          offp_data(n) = onp[i]%data(offset+k)
        end do
      end associate
    end do
    sync all ! avoid buffer being deallocated at return while still needed by another image
#ifdef __GFORTRAN__
    !deallocate(onp)
#endif
#undef __DATA_TYPE__