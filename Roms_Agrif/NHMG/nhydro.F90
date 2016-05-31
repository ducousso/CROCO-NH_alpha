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

    call define_matrices(dx, dy, h)

  end subroutine nhydro_init

  !--------------------------------------------------------------
  subroutine nhydro_solve(u,v,w)

    real(kind=rp), dimension(:,:,:), pointer, intent(inout) :: u,v,w

    real(kind=rp)    :: tol = 1.e-12
    integer(kind=ip) :: maxite = 50

    grid(1)%p(:,:,:) = 0._rp

    call solve(tol,maxite)

    if (netcdf_output) then
       call write_netcdf(grid(1)%p,vname='p',netcdf_file_name='p_end.nc',rank=myrank)
       call write_netcdf(grid(1)%r,vname='r',netcdf_file_name='r_end.nc',rank=myrank)
    endif

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
