    class(index_map), intent(in) :: this
    __DATA_TYPE__, intent(in), target :: src(:,:,:)
    __DATA_TYPE__, intent(inout) :: dest(:,:,:)
        
    type box
      __DATA_TYPE__, pointer :: data(:,:,:)
    end type
    type(box), allocatable :: buffer[:]

    ASSERT(size(src,3) >= merge(this%global_size,0,this%is_root))
    ASSERT(size(dest,3) >= this%onp_size)
    
    allocate(buffer[*])
    if (this_image() == this%root) buffer%data => src
    sync all
    dest(:,:,:this%onp_size) = buffer[this%root]%data(:,:,this%first_gid:this%last_gid)
    sync all
#undef __DATA_TYPE__
