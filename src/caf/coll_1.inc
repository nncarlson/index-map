    class(index_map), intent(in) :: this
    __DATA_TYPE__, intent(in) :: src(:)
    __DATA_TYPE__, intent(inout), target :: dest(:)
        
    type box
      __DATA_TYPE__, pointer :: data(:)
    end type
    type(box), allocatable :: buffer[:]

    ASSERT(size(dest) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(src) >= this%onp_size)
    
    allocate(buffer[*])
    if (this_image() == this%root) buffer%data => dest
    sync all
    buffer[this%root]%data(this%first_gid:this%last_gid) = src(:this%onp_size)
    sync all
#undef __DATA_TYPE__
