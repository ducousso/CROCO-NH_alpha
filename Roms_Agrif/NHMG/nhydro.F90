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
    real(kind=rp), dimension(1:nx+1,0:ny+1,1:nz),   intent(inout) :: u
    real(kind=rp), dimension(0:nx+1,1:ny+1,1:nz),   intent(inout) :: v
    real(kind=rp), dimension(0:nx+1,0:ny+1,0:nz), intent(inout) :: w

    real(kind=rp)    :: tol = 1.e-12
    integer(kind=ip) :: maxite = 10

    real(kind=rp), dimension(:,:),   pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zr,zw
    real(kind=rp), dimension(:,:,:), pointer :: dz,dzw
    real(kind=rp), dimension(:,:,:), pointer :: Arx, Ary
    real(kind=rp), dimension(:,:),   pointer :: Arz 
    real(kind=rp), dimension(:,:,:), pointer :: zy,zx,zyw,zxw,zydx,zxdy
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

    write(*,*) 'l u', lbound(u), 'u u', ubound(u)
    write(*,*) 'l v', lbound(v), 'u v', ubound(v)
    write(*,*) 'l w', lbound(w), 'u w', ubound(w)

    nh = grid(1)%nh

    dx => grid(1)%dx
    dy => grid(1)%dy
    zr => grid(1)%zr
    zw => grid(1)%zw

    !! Cell heights
    allocate(dz(nz,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             dz(k,j,i) = zw(k+1,j,i)-zw(k,j,i) !!  cell height at rho-points
          enddo
       enddo
    enddo
    allocate(dzw(nz+1,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          dzw(1,j,i) = zr(1,j,i)-zw(1,j,i) !!
          do k = 2,nz
             dzw(k,j,i) = zr(k,j,i)-zr(k-1,j,i) !!  cell height at w-points
          enddo
          dzw(nz+1,j,i) = zw(nz+1,j,i)-zr(nz,j,i) !!
       enddo
    enddo

    !!  Areas
    allocate(Arx(nz,ny,nx+1))
    do i = 1,nx+1
       do j = 1,ny
          do k = 1,nz
             Arx(k,j,i) = 0.25_8*(dz(k,j,i)+dz(k,j,i-1))*(dy(j,i)+dy(j,i-1)) 
          enddo
       enddo
    enddo
    allocate(Ary(nz,ny+1,nx))
    do i = 1,nx
       do j = 1,ny+1
          do k = 1,nz
             Ary(k,j,i) = 0.25_8*(dz(k,j,i)+dz(k,j-1,i))*(dx(j,i)+dx(j-1,i)) 
          enddo
       enddo
    enddo
    allocate(Arz(0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          Arz(j,i) = dx(j,i)*dy(j,i)
       enddo
    enddo

    !! Slopes in x- and y-direction defined at rho-points
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

    allocate(zxdy(nz,0:ny+1,0:nx+1))
    allocate(zydx(nz,0:ny+1,0:nx+1))
    do k = 1,nz
       zydx(k,:,:) = zy(k,:,:)*dx(:,:)
       zxdy(k,:,:) = zx(k,:,:)*dy(:,:)
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

          k = 1 !lower level
!!$             uf(k,j,i) = u(i,j,k) &
!!$                  - 0.25*( zx(k,j,i  )*w(i,j,k-1  ) + zx(k,j,i  )*w(i,j,k  ) &
!!$                         + zx(k,j,i-1)*w(i-1,j,k-1) + zx(k,j,i-1)*w(i-1,j,k) )
          uf(k,j,i) = Arx(k,j,i)*u(i,j,k) &
               - 0.25*( zxdy(k,j,i  )*dzw(k,j,i)*w(i,j,k-1  ) + zxdy(k,j,i  )*dzw(k+1,j,i)*w(i,j,k  ) &
                      + zxdy(k,j,i-1)*dzw(k,j,i-1)*w(i-1,j,k-1) + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(i-1,j,k) )

!!$             vf(k,j,i) = v(i,j,k) &
!!$                  - 0.25*( zy(k,j  ,i)*w(i,j,k-1) + zy(k,j  ,i)*w(i,j,k) &
!!$                         + zy(k,j-1,i)*w(i,j-1,k-1) + zy(k,j-1,i)*w(i,j-1,k) )
          vf(k,j,i) =  Ary(k,j,i)*v(i,j,k) &
               - 0.25*( zydx(k,j  ,i)*dzw(k,j,i)*w(i,j,k-1) + zydx(k,j  ,i)*dzw(k+1,j,i)*w(i,j,k) &
                      + zydx(k,j-1,i)*dzw(k,j-1,i)*w(i,j-1,k-1) + zydx(k,j-1,i)*dzw(k+1,j-1,i)*w(i,j-1,k) )

          wf(k,j,i) = 0._8

          do k = 2,nz !interior levels
!!$             uf(k,j,i) = u(i,j,k) &
!!$                  - 0.25*( zx(k,j,i  )*w(i,j,k-1  ) + zx(k,j,i  )*w(i,j,k  ) &
!!$                         + zx(k,j,i-1)*w(i-1,j,k-1) + zx(k,j,i-1)*w(i-1,j,k) )
          uf(k,j,i) = Arx(k,j,i)*u(i,j,k) &
               - 0.25*( zxdy(k,j,i  )*dzw(k,j,i)*w(i,j,k-1  ) + zxdy(k,j,i  )*dzw(k+1,j,i)*w(i,j,k  ) &
                      + zxdy(k,j,i-1)*dzw(k,j,i-1)*w(i-1,j,k-1) + zxdy(k,j,i-1)*dzw(k+1,j,i-1)*w(i-1,j,k) )

!!$             vf(k,j,i) = v(i,j,k) &
!!$                  - 0.25*( zy(k,j  ,i)*w(i,j,k-1) + zy(k,j  ,i)*w(i,j,k) &
!!$                         + zy(k,j-1,i)*w(i,j-1,k-1) + zy(k,j-1,i)*w(i,j-1,k) )
          vf(k,j,i) =  Ary(k,j,i)*v(i,j,k) &
               - 0.25*( zydx(k,j  ,i)*dzw(k,j,i)*w(i,j,k-1) + zydx(k,j  ,i)*dzw(k+1,j,i)*w(i,j,k) &
                      + zydx(k,j-1,i)*dzw(k,j-1,i)*w(i,j-1,k-1) + zydx(k,j-1,i)*dzw(k+1,j-1,i)*w(i,j-1,k) )

!!$             wf(k,j,i) = w(i,j,k-1)*(1._8 + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i)) &
!!$                  - 0.25*( zx(k,j,i  )*u(i,j,k  ) + zx(k,j,i  )*u(i+1,j,k ) &
!!$                                + zx(k-1,j,i)*u(i,j,k-1) + zx(k-1,j,i)*u(i+1,j,k-1) ) &
!!$                  - 0.25*( zy(k,j,i  )*v(i,j,k  ) + zy(k,j  ,i)*v(i,j+1,k) &
!!$                         + zy(k-1,j,i)*v(i,j,k-1) + zy(k-1,j,i)*v(i,j+1,k-1) )                 
             wf(k,j,i) = Arz(j,i)*w(i,j,k-1)*(1._8 + zxw(k,j,i)*zxw(k,j,i)+zyw(k,j,i)*zyw(k,j,i)) &
                  - 0.25*( zxdy(k,j,i  )*dx(j,i)*u(i,j,k  ) + zxdy(k,j,i  )*dx(j,i+1)*u(i+1,j,k ) &
                         + zxdy(k-1,j,i)*dx(j,i)*u(i,j,k-1) + zxdy(k-1,j,i)*dx(j,i+1)*u(i+1,j,k-1) ) &
                  - 0.25*( zydx(k,j,i  )*dy(j,i)*v(i,j,k  ) + zydx(k,j  ,i)*dy(j+1,i)*v(i,j+1,k) &
                         + zydx(k-1,j,i)*dy(j,i)*v(i,j,k-1) + zydx(k-1,j,i)*dy(j+1,i)*v(i,j+1,k-1) )
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
       !!
       !!  um = um - dxp
       !!  
       !!  
       !!  care : - arrays shape

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
