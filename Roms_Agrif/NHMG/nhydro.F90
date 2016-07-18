#include "cppdefs.h"
#ifdef NHMG
module nhydro

  use mg_mpi
  use mg_grids
  use mg_namelist
  use mg_tictoc
  use mg_mpi_exchange
  use mg_netcdf_out
  use mg_compute_rhs
  use mg_correct_uvw
  use mg_solvers
  use mg_compute_barofrc


  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhydro_init(  &
       nx, ny, nz,         &
       npxg, npyg,         &
       dx,dy,h,            &
       hc,theta_b,theta_s, &
       test)

    integer(kind=ip), intent(in) :: nx, ny, nz
    integer(kind=ip), intent(in) :: npxg, npyg
    real(kind=rp), dimension(:,:), intent(in) :: dx, dy, h
    real(kind=rp),  intent(in) :: hc, theta_b, theta_s
    character(len=*), optional :: test

    call tic(1,'nhydro_init')

    call mg_mpi_init()

    hlim      = hc
    nhtheta_b = theta_b
    nhtheta_s = theta_s

    !- read the NonHydro namelist file if it is present 
    !- else default values and print them (or not).
    call read_nhnamelist(vbrank=myrank)

    call define_grids(npxg, npyg, nx, ny, nz)

    call define_neighbours()

    call print_grids()

    if (present(test)) then
       if (trim(test)=='seamount') then
          bench = 'seamount'
       endif
    endif

    call define_matrices(dx, dy, h)

    call toc(1,'nhydro_init')

  end subroutine nhydro_init

  !--------------------------------------------------------------
  subroutine nhydro_solve(nx,ny,nz,ua,va,wa,rua,rva)

    integer(kind=ip), intent(in) :: nx, ny, nz

    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz), target, intent(inout) :: ua
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz), target, intent(inout) :: va
    real(kind=rp), dimension(0:nx+1,0:ny+1,0:nz), target, intent(inout) :: wa
    real(kind=rp), dimension(1:nx+1,0:ny+1),      target, intent(out)   :: rua
    real(kind=rp), dimension(0:nx+1,1:ny+1),      target, intent(out)   :: rva

    real(kind=rp), dimension(:,:,:), pointer :: u, v, w
    real(kind=rp), dimension(:,:)  , pointer :: ru, rv

    real(kind=rp), dimension(:,:,:), allocatable, target :: ub, vb, wb


    real(kind=rp)    :: tol
    integer(kind=ip) :: maxite
    integer(kind=ip) :: i,j,k

    call tic(1,'nhydro_solve')

    tol    = solver_prec    ! solver_prec    is defined in the namelist file
    maxite = solver_maxiter ! solver_maxiter is defined in the namelist file

    ! Reshape u,v and w array indexing
    allocate(ub(1:nz,0:ny+1,nx+1))
    allocate(vb(1:nz,ny+1,0:nx+1))
    allocate(wb(0:nz+1,0:ny+1,0:nx+1))

    do i = 1, nx+1
      do j = 0,ny+1
        do k = 1,nz
          ub(k,j,i) = ua(i,j,k)
        enddo
      enddo
    enddo

    do i = 0, nx+1
      do j = 1,ny+1
        do k = 1,nz
          vb(k,j,i) = va(i,j,k)
        enddo
      enddo
    enddo

    do i = 0, nx+1
      do j = 0,ny+1
        do k = 0,nz
          wb(k,j,i) = wa(i,j,k)
        enddo
      enddo
    enddo

    u => ub
    v => vb
    w => wb

!    u => ua
!    v => va
!    w => wa
    ru => rua
    rv => rva

    if (myrank==0) write(*,*)' nhydro_solve:'

!!$    if (myrank==0) write(*,*) 'lbound(u) ',lbound(u), 'ubound(u) ',ubound(u)
!!$    if (myrank==0) write(*,*) 'lbound(v) ',lbound(v), 'ubound(v) ',ubound(v)
!!$    if (myrank==0) write(*,*) 'lbound(w) ',lbound(w), 'ubound(w) ',ubound(w)

    !- step 1 - 
    call compute_rhs(u,v,w)

    if (netcdf_output) then
       call write_netcdf(grid(1)%b,vname='b',netcdf_file_name='b.nc',rank=myrank)
    endif

    !- step 2 -
    call solve_p(tol,maxite)

    if (netcdf_output) then
       call write_netcdf(grid(1)%p,vname='p',netcdf_file_name='p.nc',rank=myrank)
       call write_netcdf(grid(1)%r,vname='r',netcdf_file_name='r.nc',rank=myrank)
    endif

    !- step 3 -
    call correct_uvw(u,v,w)

    !- step 4 -
    call compute_barofrc(ru,rv)

   if (myrank==0) write(*,*)' kji -> ijk'
   do i = 1, nx+1
      do j = 0,ny+1
        do k = 1,nz
          ua(i,j,k) = ub(k,j,i)
        enddo
      enddo
    enddo

    do i = 0, nx+1
      do j = 1,ny+1
        do k = 1,nz
          va(i,j,k) = vb(k,j,i)
        enddo
      enddo
    enddo

    do i = 0, nx+1
      do j = 0,ny+1
        do k = 0,nz
          wa(i,j,k) = wb(k,j,i)
        enddo
      enddo
    enddo

    if (myrank==0) write(*,*)'deallocate ub, vb and wb' 
    u => null()
    v => null()
    w => null()
    deallocate(ub)
    deallocate(vb)
    deallocate(wb)
    call toc(1,'nhydro_solve')	

    if (myrank==0) write(*,*)' nhydro_solve end !!!'
 
  end subroutine nhydro_solve

  !--------------------------------------------------------------
  subroutine nhydro_clean()

    real(kind=rp) :: tstart,tend,perf

    call cpu_time(tstart)

    call grids_dealloc()

    call print_tictoc()

    call cpu_time(tend)

    if (myrank == 0) write(*,*)'nhydro_clean time:',tend-tstart

  end subroutine nhydro_clean

end module nhydro
#else
        module nhydro_empty
        end module
#endif
