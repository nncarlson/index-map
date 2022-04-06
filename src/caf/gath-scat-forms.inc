

    !! GATHER_OFFP: remote read of scattered on-process data to off-process data

    onp%data => onp_data

    n = 0
    do i = 1, num_images()
      offset = this%onp_displs(i)
      do j = 1, this%offp_counts(i)
        n = n + 1
        offp_data(n) = onp[i]%data(this%onp_index(offset+j)) !! WRONG: ONP_INDEX from image i too!
      end do
    end do

    !! GATHER_OFFP: remote read of buffered on-process data to off-process data

    do j = 1, size(this%onp_index)
      buffer%data(j) = onp_data(this%onp_index(j))
    end do

    n = 0
    do i = 1, num_images()
      offset = this%onp_displs(i)
      do j = 1, this%offp_counts(i)
        n = n + 1
        offp_data(n) = buffer[i]%data(offset+j)
      end do
    end do

    !! GATHER_OFFP: remote write to off-process data from scattered on-process data

    offp%data => offp_data

    n = 0
    do i = 1, num_images()
      offset = this%offp_displs(i)
      do j = 1, this%onp_counts(i)
        n = n + 1
        offp[i]%data(offset+j) = onp_data(this%onp_index(n))
      end do
    end do

    !! GATHER_OFFP: remote write to off-process data from buffered on-process data

    offp%data => offp_data
    onp_buffer = onp_data(onp_index)

    n = 0
    do i = 1, num_images()
      offset = this%offp_displs(i)
      do j = 1, this%onp_counts(i)
        n = n + 1
        offp[i]%data(offset+j) = onp_buffer(n)
      end do
    end do

    !! SCATTER_OFFP: remote read of off-process data to scattered on-process data

    offp%data => offp_data

    n = 0
    do i = 1, num_images()
      offset = this%offp_displs(i)
      do j = 1, this%onp_counts(i)
        n = n + 1
        onp_data(this%onp_index(n)) += offp[i]%data(offset+j)
      end do
    end do

    !! SCATTER_OFFP: remote read to buffered data

    offp%data => offp_data

    n = 0
    do i = 1, num_images()
      offset = this%offp_displs(i)
      do j = 1, this%onp_counts(i)
        n = n + 1
        onp_buffer(n) = offp[i]%data(offset+j)
      end do
    end do
    
    do j = 1, size(this%onp_index)
      onp_data(this%onp_index(j)) += onp_buffer(j)
    end do

    !! SCATTER_OFFP: remote write of off-process data to scattered on-process data

       ! This organization is not possible without some form of memory locking
       ! because multiple images may attempt to update the same on-process element
       ! at the same time. This is due to the mapping from off-process elements to
       ! corresponding on-process elements being many-to-one.

    !! SCATTER_OFFP: remote write of off-process data to buffered on-process data

    n = 0
    do i = 1, num_images()
      offset = this%onp_displs(i)
      do j = 1, this%offp_counts(i)
        n = n + 1
        onp_buffer[i]%data(offset+j) = offp_data(n)
      end do
    end do
    
    do j = 1, size(this%onp_index)
      onp_data(this%onp_index(j)) += onp_buffer%data(j)
    end do
