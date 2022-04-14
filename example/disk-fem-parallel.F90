!! DISK_FEM_PARALLEL
!!
!! This example program solves the heat equation \(u_t = \Delta u\) on the
!! unit disk subject to 0 boundary conditions. It uses a Galerkin finite
!! element discretization with bilinear elements on a regular 2D Cartesian
!! mesh that contains the unit disk. Only those cells with at least one node
!! contained inside the unit disk are included in the problem. The unknowns
!! are nodal values of \(u\). In order to avoid the need for a linear solver,
!! a lumped mass matrix and simple first-order forward Euler time stepping
!! is used to solve from a uniform initial condition to a final time. Domain
!! decomposition is used to parallelize the computation. The cells are numbered
!! according to the usual order of the elements of a rank-2 array (but only for
!! those included in the problem) and then partitioned into approximately equal
!! blocks, one per process (i.e., MPI rank or coarray image). For reference,
!! see the serial implementation disk-fem-serial.F90.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright 2022 Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!! Permission is hereby granted, free of charge, to any person obtaining a
!! copy of this software and associated documentation files (the "Software"),
!! to deal in the Software without restriction, including without limitation
!! the rights to use, copy, modify, merge, publish, distribute, sublicense,
!! and/or sell copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following conditions:
!!
!! The above copyright notice and this permission notice shall be included
!! in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!! DEALINGS IN THE SOFTWARE.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program disk_fem_parallel

  use,intrinsic ::iso_fortran_env, only: i8 => int64, r8 => real64
#ifndef USE_CAF
  use mpi
#endif
  use index_map_type
  implicit none

  integer, parameter :: NZ = 257  ! number of zones in each dimension
  integer :: ierr, nproc, rank, n, nnode, ncell, j, step, nstep
  integer, allocatable :: cell_mask(:,:), node_mask(:,:), cnode(:,:), cnode_local(:,:)
  integer, allocatable :: cell_bsizes(:), node_bsizes(:)
  real(r8), allocatable :: u(:), u_local(:), f(:)
  type(index_map) :: cell_map, node_map
  real(r8) :: dt, dx, c, tfinal, smat(4,4)
  integer(i8) :: t1, t2, rate, s(2), s1, s2

#ifdef USE_CAF
  nproc = num_images()
  rank = this_image() - 1
  if (rank == 0) write(*,'(a,i0,a)') 'Running with ', nproc, ' CAF images'
#else
  call MPI_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (rank == 0) write(*,'(a,i0,a)') 'Running with ', nproc, ' MPI ranks'
#endif

  !! The NODE_MASK and CELL_MASK arrays marks grid node and cells inside the
  !! unit disk with a unique positive index. Other nodes and cells are marked
  !! with a 0. A cell is considered inside if any of its 4 nodes is inside.
  !! NNODE and NCELL are the number of nodes and cells contained in the unit
  !! disk and included in the problem. NNODE is also the number of
  !! node-centered unknowns.
  if (rank == 0) then
    allocate(node_mask(0:NZ,0:NZ), cell_mask(NZ,NZ))
    call unit_disk_mask(node_mask, nnode, cell_mask, ncell)
  end if
#ifdef USE_CAF
  call co_broadcast(ncell, source_image=1)
  call co_broadcast(nnode, source_image=1)
#else
  call MPI_Bcast(ncell, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(nnode, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

  !! Create the indirect indexing array CNODE. CNODE(:,j) are the indices of
  !! the 4 nodes adjacent to cell j, or 0 if the node does not lie inside the
  !! unit disk.
  n = merge(ncell, 0, rank==0)
  allocate(cnode(4,n))
  if (rank == 0) call fill_cnode(cell_mask, node_mask, cnode)

  !! Create a mapping of the cell index set onto NPROC processes with
  !! nearly equal block sizes. This index set will have no overlap.
  !! Also create a mapping of the node index set onto NPROC processes.
  n = merge(nproc, 0, rank==0)
  allocate(cell_bsizes(n), node_bsizes(n))
  if (rank == 0) then
    n = ncell/nproc
    cell_bsizes = n
    n = ncell - n*nproc
    cell_bsizes(1:n) = cell_bsizes(1:n) + 1
    call partition_nodes(cnode, cell_bsizes, nnode, node_bsizes)
  end if
  call cell_map%init(cell_bsizes)
  call node_map%init(node_bsizes)

  !! Localize the CNODE indirect indexing array. This adds the overlap to the
  !! node index map that is required for closure of the localized CNODE array.
  call cell_map%localize_index_array(cnode, node_map, cnode_local)

  !! U is the global array of node-centered unknowns, and U_LOCAL is the
  !! corresponding distributed array. For convenience, an extra element
  !! is prepended to the arrays for the artificial index 0 to serve as a
  !! place to hold the boundary value.
  n = merge(nnode, 0, rank==0)
  allocate(u(0:n), u_local(0:node_map%local_size))
  u = 1  ! initial conditions
  call node_map%distribute(u(1:), u_local(1:))
  u_local(0) = 0 ! boundary value

  !! Time step out to time TFINAL using the explicit first order Euler
  !! method with a fixed stable time step DT.
  dx = 1.0_r8/NZ
  dt = 0.25_r8*dx**2
  c  = dt/(6*dx**2)

  !! FE stiffness matrix
  smat = reshape([4,-1,-2,-1, -1,4,-1,-2, -2,-1,4,-1, -1,-2,-1,4], shape=[4,4])

  tfinal = 0.05_r8
  nstep = ceiling(tfinal/dt)
  tfinal = nstep*dt

  allocate(f, mold=u_local)
  call system_clock(t1)
  s = 0
  do step = 1, nstep
    call system_clock(s1)
    call node_map%gather_offp(u_local(1:))
    call system_clock(s2)
    s(1) = s(1) + (s2-s1)
    f = 0
    do j = 1, cell_map%onp_size ! the usual FE assembly process
#ifdef COMPACT_UPDATE
      f(cnode_local(:,j)) = f(cnode_local(:,j)) + matmul(smat, u_local(cnode_local(:,j)))
#else
      !f(cnode_local(:,j)) = f(cnode_local(:,j)) + matmul(smat, u_local(cnode_local(:,j)))
      block
        real(r8) :: ucell(4), fcell(4)
        integer :: i
        do i = 1, 4
          ucell(i) = u_local(cnode_local(i,j))
        end do
        fcell = matmul(smat, ucell)
        do i = 1, 4
          f(cnode_local(i,j)) = f(cnode_local(i,j)) + fcell(i)
        end do
      end block
#endif
    end do
    call system_clock(s1)
    call node_map%scatter_offp_sum(f(1:)) ! complete the FE assembly
    call system_clock(s2)
    s(2) = s(2) + (s2-s1)
    do j = 1, node_map%onp_size
      u_local(j) = u_local(j) - c*f(j)
    end do
  end do
  call system_clock(t2, rate)

  !! Gather the distributed solution onto rank 0 and output.
  call node_map%collate(u_local(1:), u(1:))
  if (rank == 0) then
    write(*,'(a,es10.4,a)') 'Solution at t=', tfinal, ' written to out.vtk; visualize with paraview.'
    u(0) = 0 ! boundary value
    call vtk_plot(node_mask, u)
    write(*,'(g0,a,g0,1x,g0,a,2(i0,a))') (10**6)*real(t2-t1)/real(rate)/nstep, &
        ' µsec/time step (', (10**6)*real(s)/real(rate)/nstep, ' comm); ', &
        ncell/nproc, ' cells/process, ', nstep, ' steps'
  end if

#ifndef USE_CAF
  call MPI_Finalize(ierr)
#endif

contains

  !! Generate a grid mask arrays for a unit disk. A node_mask element for a
  !! grid node that does not lie in the interior of the unit disk is set to
  !! zero. Other elements are assigned sequential integers starting with 1
  !! according to the normal traveral order of a rank-2 array. The number of
  !! grid nodes so contained inside the unit disk is returned in NNODE. A
  !! A cell_mask array element is set to zero if all the nodes surrounding
  !! the corresponding cell have zero node_mask values. Otherwise the elements
  !! are assigned sequential integers starting with 1 according to the normal
  !! traversal order of a rank-2 array. The number of grid cells so numbered
  !! is returned in NCELL.

  subroutine unit_disk_mask(node_mask, nnode, cell_mask, ncell)

    integer, intent(out) :: node_mask(0:,0:), nnode, cell_mask(:,:), ncell

    real :: dx, dy, x, y
    integer :: n, i, j

    dx = 2.0/NZ
    dy = 2.0/NZ
    n = 0
    y = -1
    do j = 0, NZ
      if (j == NZ) y = 1
      x = -1
      do i = 0, NZ
        node_mask(i,j) = 0
        if (i == NZ) x = 1
        if (x**2 + y**2 < 1) then
          n = n + 1
          node_mask(i,j) = n
        end if
        x = x + dx
      end do
      y = y + dy
    end do
    nnode = n

    n = 0
    do j = 1, NZ
      do i = 1, NZ
        cell_mask(i,j) = 0
        if (node_mask(i-1,j-1) == 0) then
          if (node_mask(i,j-1) == 0) then
            if (node_mask(i-1,j) == 0) then
              if (node_mask(i,j) == 0) cycle
            end if
          end if
        end if
        n = n + 1
        cell_mask(i,j) = n
      end do
    end do
    ncell = n

  end subroutine

  !! Fill the cell node indirect indexing array using the node mask array to
  !! obtain the adjacent node indices (or 0 where the node is not included).

  subroutine fill_cnode(cell_mask, node_mask, cnode)
    integer, intent(in) :: cell_mask(:,:), node_mask(0:,0:)
    integer, intent(out) :: cnode(:,:)
    integer :: n, i, j
    n = 0
    do j = 1, NZ
      do i = 1, NZ
        n = cell_mask(i,j)
        if (n > 0) then
          cnode(1,n) = node_mask(i-1,j-1)
          cnode(2,n) = node_mask(i,j-1)
          cnode(3,n) = node_mask(i,j)
          cnode(4,n) = node_mask(i-1,j)
        end if
      end do
    end do
  end subroutine

  !! Generate a block partitioning of the nodes that conforms to the given
  !! block partitioning of the cells: a node is assigned to a partition of
  !! one of its adjacent cells. For nodes with multiple choices, a simple
  !! one is to choose the partition with the lowest index (greedy). In
  !! general, the nodes may need to be renumbered in order to for this to
  !! produce the required block partitioning. However in this specific case
  !! and for the way cells and nodes are numbered, no renumbering appears
  !! to be required. This is checked, and the program will stop if it is
  !! found to be otherwise.

  subroutine partition_nodes(cnode, cell_bsizes, num_node, node_bsizes)
    integer, intent(in) :: cnode(:,:), cell_bsizes(:), num_node
    integer, intent(out) :: node_bsizes(:)
    integer :: part(num_node), i, j, k, n, i1
    part = 0
    n = 0
    do j = 1, size(cell_bsizes)
      do i = 1, cell_bsizes(j)
        n = n + 1
        do k = 1, size(cnode,1)
          if (cnode(k,n) == 0) cycle
          if (part(cnode(k,n)) == 0) part(cnode(k,n)) = j
        end do
      end do
    end do
    if (minval(part) == 0) error stop 1
    if (any(part(2:) < part(1:size(part)-1))) error stop 2
    i1 = 1
    do j = 1, size(node_bsizes)
      do i = i1, size(part)
        if (part(i) > j) exit
      end do
      node_bsizes(j) = i - i1
      i1 = i
    end do
  end subroutine

  !! Write a VTK plot file of the solution U. The solution is embedded in a
  !! rectangular grid with boundary values (0) for grid points outside the
  !! masked region.

  subroutine vtk_plot(mask, u)
    integer, intent(in) :: mask(0:,0:)
    real(r8), intent(in) :: u(0:)
    integer :: lun, i, j
    real, allocatable :: u_plot(:,:)
    open(newunit=lun,file='out.vtk')
    write(lun,'("# vtk DataFile Version 3.0")')
    write(lun,'("example solution")')
    write(lun,'("ASCII")')
    write(lun,'("DATASET STRUCTURED_POINTS")')
    write(lun,'("DIMENSIONS",3(1x,i0))') NZ+1, NZ+1, 1
    write(lun,'("ORIGIN 0 0 0")')
    write(lun,'("SPACING",3(1x,g0))') 1.0/NZ, 1.0/NZ, 1
    write(lun,'("POINT_DATA ",i0)') (NZ+1)**2
    write(lun,'("SCALARS u float 1")')
    write(lun,'("LOOKUP_TABLE default")')
    allocate(u_plot(0:NZ,0:NZ))
    do j = 0, NZ
      do i = 0, NZ
        u_plot(i,j) = u(mask(i,j))
      end do
    end do
    write(lun,'(g0)') u_plot
    close(lun)
  end subroutine

end program
