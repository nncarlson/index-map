!! DISK_FEM_SERIAL
!!
!! This example program solves the heat equation \(u_t = \Delta u\) on the
!! unit disk subject to 0 boundary conditions. It uses a Galerkin finite
!! element discretization with bilinear elements on a regular 2D Cartesian
!! mesh that contains the unit disk. Only those cells with at least one node
!! contained inside the unit disk are included in the problem. The unknowns
!! are nodal values of \(u\). In order to avoid the need for a linear solver,
!! we use a lumped mass matrix and use simple first-order forward Euler time
!! stepping to solve from a uniform initial condition to a final time. This
!! is a serial implementation intended as a reference for the parallel
!! implementation disk-fem-parallel.F90.
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

program disk_fem_serial

  use,intrinsic ::iso_fortran_env, only: i8 => int64, r8 => real64
  implicit none

  integer, parameter :: NZ = 257  ! number of zones in each dimension
  integer :: nnode, ncell, j, step, nstep
  integer, allocatable :: cell_mask(:,:), node_mask(:,:), cnode(:,:)
  real(r8), allocatable :: u(:), f(:)
  real(r8) :: dt, dx, c, tfinal, S(4,4)
  integer(i8) :: t1, t2, rate

  !! The NODE_MASK and CELL_MASK arrays marks grid node and cells inside the
  !! unit disk with a unique positive index. Other nodes and cells are marked
  !! with a 0. A cell is considered inside if any of its 4 nodes is inside.
  !! NNODE and NCELL are the number of nodes and cells contained in the unit
  !! disk and included in the problem. NNODE is also the number of
  !! node-centered unknowns.
  allocate(node_mask(0:NZ,0:NZ), cell_mask(NZ,NZ))
  call unit_disk_mask(node_mask, nnode, cell_mask, ncell)

  !! Create the indirect indexing array CNODE. CNODE(:,j) are the indices of
  !! the 4 nodes adjacent to cell j, or 0 if the node does not lie inside the
  !! unit disk.
  allocate(cnode(4,ncell))
  call fill_cnode(cell_mask, node_mask, cnode)

  !! U is the array of node-centered unknowns. For convenience, an extra element
  !! is prepended to the array for the artificial index 0 to serve as a place to
  !! hold the boundary value.
  allocate(u(0:nnode))
  u = 1  ! initial conditions
  u(0) = 0 ! boundary value

  !! Time step out to time TFINAL using the explicit first order Euler
  !! method with a fixed stable time step DT.
  dx = 1.0_r8/NZ
  dt = 0.25_r8*dx**2
  c  = dt/(6*dx**2)

  !! FE stiffness matrix
  S = reshape([4,-1,-2,-1, -1,4,-1,-2, -2,-1,4,-1, -1,-2,-1,4], shape=[4,4])

  tfinal = 0.05_r8
  nstep = ceiling(tfinal/dt)
  tfinal = nstep*dt

  allocate(f, mold=u)
  call system_clock(t1)
  do step = 1, nstep
    f = 0
    do j = 1, ncell ! the usual FE assembly process
#ifdef COMPACT_UPDATE
      f(cnode(:,j)) = f(cnode(:,j)) + matmul(S, u(cnode(:,j)))
#else
      !f(cnode(:,j)) = f(cnode(:,j)) + matmul(S, u(cnode(:,j)))
      block
        real(r8) :: ucell(4), z, fcell(4)
        integer :: i
        do i = 1, 4
          ucell(i) = u(cnode(i,j))
        end do
        fcell = matmul(S,ucell)
        do i = 1, 4
          f(cnode(i,j)) = f(cnode(i,j)) + fcell(i)
        end do
      end block
#endif
    end do
    f(0) = 0
    u = u - c*f
  end do
  call system_clock(t2, rate)

  !! Output the solution.
  write(*,'(a,es10.4,a)') 'Solution at t=', tfinal, ' written to out.vtk; visualize with paraview.'
  call vtk_plot(node_mask, u)
  write(*,'(g0,a,2(i0,a))') (10**6)*real(t2-t1)/real(rate)/nstep, &
      ' µsec per time step; ', ncell, ' cells, ', nstep, ' steps'

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
