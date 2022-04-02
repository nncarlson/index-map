!! DISK_FV_SERIAL
!!
!! This example program solves the heat equation \(u_t = \Delta u\) on the
!! unit disk subject to 0 boundary conditions. It uses a simple finite volume
!! discretization on a regular 2D Cartesian mesh that contains the unit disk.
!! Only those cells whose center is contained in the unit disk are included
!! in the problem. The unknowns are cell-centered values of \(u\). Simple
!! first-order forward Euler time stepping is used to solve from a uniform
!! initial condition to a final time. This is a serial implementation intended
!! as a reference for the MPI-parallel implementation disk-fv-parallel.F90.
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

program disk_fv_serial

  use,intrinsic ::iso_fortran_env, only: i8 => int64, r8 => real64
  implicit none

  integer, parameter :: NZ = 257  ! number of zones in each dimension
  integer :: ncell, j, step, nstep, nstep0
  integer, allocatable :: mask(:,:), cnhbr(:,:)
  real(r8), allocatable :: u(:), u_prev(:)
  real(r8) :: dt, dx, c, tfinal
  integer(i8) :: t1, t2, rate

  !! The MASK array marks grid cells inside the unit disk with a unique
  !! positive index. Other cells are marked with a 0. It is sized to include
  !! a border of cells around the actual grid. NCELL is the number of cells
  !! contained in the unit disk and the number of cell-centered unknowns.
  allocate(mask(0:NZ+1,0:NZ+1))
  call unit_disk_mask(mask, ncell)

  !! Create the indirect indexing array CNHBR. CNHBR(:,j) are the indices
  !! of the 4 cells adjacent to cell j, or 0 if no neighbor.
  allocate(cnhbr(4,ncell))
  call fill_cnhbr(mask, cnhbr)

  !! U is the array of cell-centered unknowns. For convenience, an extra element
  !! is prepended to the array for the artificial index 0 to serve as a place to
  !! hold the boundary value.
  allocate(u(0:ncell))
  u = 1  ! initial conditions
  u(0) = 0 ! boundary value

  !! Time step out to time TFINAL using the explicit first order Euler
  !! method with a fixed stable time step DT.
  dx = 1.0_r8/NZ
  dt = 0.25_r8*dx**2
  c  = dt/dx**2

  tfinal = 0.05_r8
  nstep = ceiling(tfinal/dt)
  tfinal = nstep*dt
  
  nstep0 = max(1, nstep/10) ! warm-up: start timing on this step

  allocate(u_prev, mold=u)
  do step = 1, nstep
    if (step == nstep0) call system_clock(t1)
    u_prev = u
    do j = 1, ncell
      u(j) = u_prev(j) + c*(sum(u_prev(cnhbr(:,j))) - 4*u_prev(j))
    end do
  end do
  call system_clock(t2, rate)

  !! Output the solution.
  write(*,'(a,es10.4,a)') 'Solution at t=', tfinal, ' written to out.vtk; visualize with paraview.'
  call vtk_plot(mask, u)
  write(*,'(g0,a,2(i0,a))') (10**6)*real(t2-t1)/real(rate)/(nstep-nstep0+1), &
      ' µsec per time step; ', ncell, ' cells, ', nstep, ' steps'

contains

  !! Generate a grid mask array for a unit disk. A mask element corresponding
  !! to a grid cell whose center lies outside the unit disk is set to zero.
  !! Other elements are assigned sequential integers starting with 1 according
  !! to the normal traveral order of a rank-2 array. The number of grid cells
  !! so contained by the unit disk is returned in N.

  subroutine unit_disk_mask(mask, n)
    integer, intent(out) :: mask(0:,0:), n
    real :: dx, dy, x, y
    integer :: i, j
    mask = 0
    dx = 2.0/NZ
    dy = 2.0/NZ
    n = 0
    y = -1 + dy/2
    do j = 1, NZ
      x = -1 + dx/2
      do i = 1, NZ
        if (x**2 + y**2 <= 1) then
          n = n + 1
          mask(i,j) = n
        end if
        x = x + dx
      end do
      y = y + dy
    end do
  end subroutine

  !! Fill the cell neighor indirect indexing array using the mask array to
  !! obtain the cell indices of the neighbors (or 0 where there is no neighbor).

  subroutine fill_cnhbr(mask, cnhbr)
    integer, intent(in) :: mask(0:,0:)
    integer, intent(out) :: cnhbr(:,:)
    integer :: n, i, j
    n = 0
    do j = 1, NZ
      do i = 1, NZ
        n = mask(i,j)
        if (n > 0) then
          cnhbr(1,n) = mask(i,j-1)
          cnhbr(2,n) = mask(i+1,j)
          cnhbr(3,n) = mask(i,j+1)
          cnhbr(4,n) = mask(i-1,j)
        end if
      end do
    end do
  end subroutine

  !! Write a VTK plot file of the solution U. The solution is embedded in a
  !! rectangular grid with boundary values (0) for grid cells outside the
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
    write(lun,'("CELL_DATA ",i0)') NZ**2
    write(lun,'("SCALARS u float 1")')
    write(lun,'("LOOKUP_TABLE default")')
    allocate(u_plot(NZ,NZ))
    do j = 1, NZ
      do i = 1, NZ
        u_plot(i,j) = u(mask(i,j))
      end do
    end do
    write(lun,'(g0)') u_plot
    close(lun)
  end subroutine

end program
