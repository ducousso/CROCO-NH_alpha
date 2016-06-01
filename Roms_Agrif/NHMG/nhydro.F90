#include "cppdefs.h"
#ifdef NHMG
module nhydro

  use mg_mpi
  use mg_grids
  use mg_namelist
  use mg_define_matrix
  use mg_solvers
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains

  !--------------------------------------------------------------
  subroutine nhydro_init(nx,ny,nz,npxg,npyg,dx,dy,h,hc,theta_b,theta_s)

    integer(kind=ip), intent(in) :: nx, ny, nz
    integer(kind=ip), intent(in) :: npxg, npyg

    real(kind=rp), dimension(:,:), intent(in) :: dx, dy, h
    real(kind=rp),  intent(in) :: hc, theta_b, theta_s

    integer(kind=ip) :: inc=1

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

    call define_matrices(dx,dy,h)

  end subroutine nhydro_init

  !--------------------------------------------------------------
  subroutine nhydro_solve(nx,ny,nz,u,v,w)

    integer(kind=ip), intent(in) :: nx, ny, nz
    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz), intent(inout) :: u
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz), intent(inout) :: v
    real(kind=rp), dimension(0:nx+1,0:ny+1,0:nz), intent(inout) :: w

    real(kind=rp)    :: tol = 1.e-12
    integer(kind=ip) :: maxite = 10

    real(kind=rp), dimension(:,:),   pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer :: zy,zx,zyw,zxw
    real(kind=rp), dimension(:,:,:), pointer :: uf,vf,wf

!!$    integer(kind=ip) :: nz,ny,nx
    integer(kind=ip) :: nh
    integer(kind=ip) :: i,j,k

    ! compute nhydro rhs
    !!
    !!  uf =  um - wm*zx
    !!  vf =  vm - wm*zy
    !!  wf = -zx*um - zy*vm + (1+zx^2+zy^2)*wm
    !!
    !!  care : - arrays are reshaped from (i,j,k) to (k,j,i)
    !!         - CROCO z is indexed from 0 to nz while we choose 
    !!           index from 1 to nz+1 in NHMG 

    write(*,*)'l u:', lbound(u), ' u u:', ubound(u)

!!$    nx = grid(1)%nx
!!$    ny = grid(1)%ny
!!$    nz = grid(1)%nz
    nh = grid(1)%nh

    dx => grid(1)%dx
    dy => grid(1)%dy
    zr => grid(1)%zr
    zw => grid(1)%zw

    allocate(zy(nz,0:ny+1,0:nx+1))
    allocate(zx(nz,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             zy(k,j,i) = 0.5_8*(zr(k,j+1,i)-zr(k,j-1,i))/dy(j,i)
             zx(k,j,i) = 0.5_8*(zr(k,j,i+1)-zr(k,j,i-1))/dx(j,i)
          enddo
       enddo
    enddo

    allocate(zyw(nz+1,0:ny+1,0:nx+1))
    allocate(zxw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz+1
             zyw(k,j,i) = 0.5_8*(zw(k,j+1,i)-zw(k,j-1,i))/dy(j,i)
             zxw(k,j,i) = 0.5_8*(zw(k,j,i+1)-zw(k,j,i-1))/dx(j,i)
          enddo
       enddo
    enddo


    allocate(uf(nz,0:ny+1,0:nx+1))
    allocate(vf(nz,0:ny+1,0:nx+1))
    allocate(wf(nz+1,0:ny+1,0:nx+1))
    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             uf(k,j,i) = u(i,j,k) &
                  - 0.25*(zx(k,j,i  )*w(i,j,k-1  ) + zx(k,j,i  )*w(i,j,k  ) &
                  + zx(k,j,i-1)*w(i-1,j,k-1) + zx(k,j,i-1)*w(i-1,j,k) )

             vf(k,j,i) = v(i,j,k) &
                  - 0.25*(zy(k,j  ,i)*w(i,j,k-1) + zy(k,j  ,i)*w(i,j,k) &
                  + zy(k,j-1,i)*w(i,j-1,k-1) + zy(k,j-1,i)*w(i,j-1,k) )

             wf(k,j,i) = - 0.25*(zx(k,j,i  )*u(i,j,k  ) + zx(k,j,i  )*u(i+1,j,k ) &
                  + zx(k-1,j,i)*u(i,j,k-1) + zx(k-1,j,i)*u(i+1,j,k-1) ) &
                  - 0.25*(zy(k,j,i  )*v(i,j,k  ) + zy(k,j  ,i)*v(i,j+1,k) &
                  + zy(k-1,j,i)*v(i,j,k-1) + zy(k-1,j,i)*v(i,j+1,k-1) ) &
                  + w(i,j,k-1)*(1._8 + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i))
          enddo
          wf(nz+1,j,i) = 0._8
       enddo
    enddo

    call fill_halo(1,uf)
    call fill_halo(1,vf)

    do i = 1,nx
       do j = 1,ny
          do k = 1,nz
             grid(1)%b(k,j,i) =  uf(k,j,i+1)-uf(k,j,i) &
                               + vf(k,j+1,i)-vf(k,j,i) &
                               + wf(k+1,j,i)-wf(k,j,i)
          enddo
       enddo
    enddo

    ! solve for nhydro pressure
    grid(1)%p(:,:,:) = 0._rp
    call solve(tol,maxite)

    if (netcdf_output) then
       call write_netcdf(grid(1)%p,vname='p',netcdf_file_name='p_end.nc',rank=myrank)
       call write_netcdf(grid(1)%r,vname='r',netcdf_file_name='r_end.nc',rank=myrank)
    endif

    ! correct u,v,w

  end subroutine nhydro_solve

  !--------------------------------------------------------------
  subroutine nhydro_clean()

    call grids_dealloc()

  end subroutine nhydro_clean

end module nhydro
#else
        module nhydro_empty
        end module
#endif
