#include "cppdefs.h"
#ifdef NHMG
module mg_compute_barofrc

  use mg_mpi
  use mg_tictoc
  use mg_namelist
  use mg_grids
  use mg_mpi_exchange
  use mg_netcdf_out

  implicit none

contains
  !-------------------------------------------------------------------------     
  subroutine compute_barofrc(ru,rv)

    real(kind=rp), dimension(:,:)  , pointer, intent(out) :: ru,rv

    integer(kind=ip):: k, j, i
    integer(kind=ip):: nx, ny, nz

    real(kind=rp), dimension(:,:)  , pointer :: dx,dy
    real(kind=rp), dimension(:,:,:), pointer :: zw
    real(kind=rp), dimension(:,:)  , pointer :: dxu,dyv
    real(kind=rp), dimension(:,:,:), pointer :: dz
    real(kind=rp), dimension(:,:,:), pointer :: p

    !NG comment: constants in a mg_cst.f90 file ?
    real(kind=rp), parameter :: hlf  = 0.5_rp

    nx = grid(1)%nx
    ny = grid(1)%ny
    nz = grid(1)%nz

    dx => grid(1)%dx
    dy => grid(1)%dy
    zw => grid(1)%zw

    if (myrank==0) write(*,*)'- compute barofrc:'

    !! Cell heights
    allocate(dz(nz,0:ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 0,ny+1
          do k = 1,nz
             dz(k,j,i) = zw(k+1,j,i)-zw(k,j,i)
          enddo
       enddo
    enddo

    !! Cell widths
    allocate(dxu(0:ny+1,nx+1))
    do i = 1,nx+1
       do j = 0,ny+1
          dxu(j,i) = hlf * (dx(j,i)+dx(j,i-1))
       enddo
    enddo
    allocate(dyv(ny+1,0:nx+1))
    do i = 0,nx+1
       do j = 1,ny+1
          dyv(j,i) = hlf * (dy(j,i)+dy(j-1,i))
       enddo
    enddo

    !! Compute
    p => grid(1)%p

    ru(:,:) = 0._8
    rv(:,:) = 0._8

    do i = 1,nx+1
       do j = 0,ny+1 
          do k = 1,nz
             ru(i,j) = ru(i,j) - hlf / dxu(j,i) *(dz(k,j,i)+dz(k,j,i-1)) *(p(k,j,i)-p(k,j,i-1))
          enddo
       enddo
    enddo

    do i = 0,nx+1
       do j = 1,ny+1 
          do k = 1,nz
             rv(i,j) = rv(i,j) - hlf / dyv(j,i) *(dz(k,j,i)+dz(k,j-1,i)) *(p(k,j,i)-p(k,j-1,i))
          enddo
       enddo
    enddo

!!$k=N
!!$ rsurf(i,j)=rho(i,j,N) + (rho(i,j,N)-rho(i,j,N-1))
!!$     &                             *(z_w(i,j,N)-z_r(i,j,N))
!!$     &                           /(z_r(i,j,N)-z_r(i,j,N-1))
!!$ pgrd(i)=(g+cff*(rsurf(i-1,j)+rsurf(i,j)))*( z_w(i-1,j,N)
!!$     &                                                   -z_w(i,j,N))
!!$ 
!!$     &     +cff*( (rho(i-1,j,N)-rsurf(i,j))*(z_w(i-1,j,N)-z_r(i,j,N))
!!$     &           +(rsurf(i-1,j)-rho(i,j,N))*(z_w(i,j,N)-z_r(i-1,j,N))
!!$     &                                                              )
!!$ ru(i,j,N) = 0.5*(HZR(i,j,N)+HZR(i-1,j,N)) *on_u(i,j) *pgrd(i)
!!$
!!$do k=N-1,1,-1
!!$ pgrd(i)=pgrd(i)-cff*(
!!$ &                                   (rho(i,j,k+1)-rho(i-1,j,k))
!!$     &                                  *(z_r(i-1,j,k+1)-z_r(i,j,k))
!!$     &                                  +(rho(i,j,k)-rho(i-1,j,k+1))
!!$     &                                  *(z_r(i,j,k+1)-z_r(i-1,j,k))
!!$     &                                                             )
!!$ ru(i,j,k) = 0.5*(HZR(i,j,k)+HZR(i-1,j,k)) *on_u(i,j) *pgrd(i)

    deallocate(dz)
    deallocate(dxu)
    deallocate(dyv)

  end subroutine compute_barofrc

end module mg_compute_barofrc
#else
        module mg_compute_barofrc_empty
        end module
#endif
