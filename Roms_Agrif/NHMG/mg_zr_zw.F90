#include "cppdefs.h"
#ifdef NHMG
module mg_zr_zw

  implicit none

  interface setup_zr_zw
     module procedure           &
          setup_zr_zw_seamount, &
          setup_zr_zw_croco
  end interface setup_zr_zw

contains

  !-------------------------------------------------------------------------     
  subroutine setup_zr_zw_seamount(h, zr, zw)

    integer(kind=4),parameter :: ip=4,rp=8
    real(kind=rp), dimension(:,:)  , pointer, intent(in)  :: h
    real(kind=rp), dimension(:,:,:), pointer, intent(out) :: zr
    real(kind=rp), dimension(:,:,:), pointer, intent(out) :: zw

    integer(kind=ip):: nx,ny,nz,nh

    integer(kind=ip):: i,j,k

    nz = size(zr,dim=1)
    nh = abs(lbound(zr,dim=2)-1) ! Assume that "one" is the first indice !
    !                            !  1 -> nh = 0  <= abs( 1-1) = 0
    !                            !  0 -> nh = 1  <= abs( 0-1) = 1
    !                            ! -1 -> nh = 2  <= abs(-1-1) = 2
    ny = size(zr,dim=2) - 2*nh
    nx = size(zr,dim=3) - 2*nh

    do i = 1-nh,nx+nh
       do j = 1-nh,ny+nh
          do k = 1,nz
             zr(k,j,i) = (real(k,kind=rp)-0.5_rp)*h(j,i)/real(nz,kind=rp) - h(j,i)
             zw(k,j,i) = (real(k,kind=rp)-1.0_rp)*h(j,i)/real(nz,kind=rp) - h(j,i)
          enddo
          zw(nz+1,j,i) = 0.0_rp
       enddo
    enddo

  end subroutine setup_zr_zw_seamount

  !-------------------------------------------------------------------------     
  subroutine setup_zr_zw_croco(hc,theta_b,theta_s,zeta,h,zr,zw,coord_type)

    ! compute zr and zw from zeta(j,i) and h(j,i)

    integer(kind=4),parameter :: ip=4,rp=8

    real(kind=rp), intent(in) :: hc
    real(kind=rp), intent(in) :: theta_b
    real(kind=rp), intent(in) :: theta_s
    real(kind=rp), dimension(:,:)  , pointer, intent(in)  :: zeta, h
    real(kind=rp), dimension(:,:,:), pointer, intent(out) :: zr, zw

    character(len=*), optional, intent(in) :: coord_type

    character(len=16) :: coord

    real(kind=rp), parameter :: nul  = 0._rp
    real(kind=rp), parameter :: one  = 1._rp
    real(kind=rp), parameter :: hlf  = 0.5_rp

    integer(kind=ip):: nx,ny,nz,nh

    integer(kind=ip):: i,j,k

    real(kind=rp) :: cff,cff1,cff2
    real(kind=rp) :: hinv
    real(kind=rp) :: sc_w,sc_r,cff_w,cff_r
    real(kind=rp) :: z_w0,z_r0
    real(kind=rp) :: cs_r,cs_w
    real(kind=rp) :: csrf, cswf

    if (present(coord_type)) then
       coord = trim(coord_type)
    else
       coord = 'new_s_coord'
    endif

    nz = size(zr,dim=1)
    nh = abs(lbound(zr,dim=2)-1) ! Assume that "one" is the first indice !
    !                            !  1 -> nh = 0  <= abs( 1-1) = 0
    !                            !  0 -> nh = 1  <= abs( 0-1) = 1
    !                            ! -1 -> nh = 2  <= abs(-1-1) = 2
    ny = size(zr,dim=2) - 2*nh
    nx = size(zr,dim=3) - 2*nh
    nh = 1

    ! vertical coordinate

    !---------------------!
    !- New s-coordinates -!
    !---------------------!
    if (trim(coord) == 'new_s_coord') then

       do i = 1-nh,nx+nh

          do j = 1-nh,ny+nh

             cff=one/real(nz,kind=rp)

             hinv = one / (h(j,i)+hc)

             do k = 1,nz

                sc_r = cff*(real(k-nz,kind=rp)-hlf)

                if (theta_s > nul) then
                   csrf=(one-cosh(theta_s*sc_r))/(cosh(theta_s)-one)
                else
                   csrf=-sc_r**2
                endif

                if (theta_b > nul) then
                   cs_r=(exp(theta_b*csrf)-one)/(one-exp(-theta_b))
                else
                   cs_r=csrf
                endif

                sc_w = cff*real(k-1-nz,kind=rp)

                if (theta_s > nul) then
                   cswf=(one-cosh(theta_s*sc_w))/(cosh(theta_s)-one)
                else
                   cswf=-sc_w**2
                endif

                if (theta_b > nul) then
                   cs_w=(exp(theta_b*cswf)-one)/(one-exp(-theta_b))
                else
                   cs_w=cswf
                endif

                cff_w = hc * sc_w
                cff_r = hc * sc_r

                z_w0 = cff_w + cs_w*h(j,i)
                z_r0 = cff_r + cs_r*h(j,i)

                zw(k,j,i) = z_w0*h(j,i)*hinv + zeta(j,i)*(1.+z_w0*hinv)
                zr(k,j,i) = z_r0*h(j,i)*hinv + zeta(j,i)*(1.+z_r0*hinv)

             enddo

             k = nz+1

             sc_w = cff*real(k-1-nz,kind=rp)

             if (theta_s > nul) then
                cswf=(one-cosh(theta_s*sc_w))/(cosh(theta_s)-one)
             else
                cswf=-sc_w**2
             endif

             if (theta_b > nul) then
                cs_w=(exp(theta_b*cswf)-one)/(one-exp(-theta_b))
             else
                cs_w=cswf
             endif

             cff_w = hc * sc_w
             z_w0 = cff_w + cs_w*h(j,i)
             zw(k,j,i) = z_w0*h(j,i)*hinv + zeta(j,i)*(1.+z_w0*hinv)
             
          enddo ! j
       enddo ! i

       !---------------------!
       !- Sigma coordinates -!
       !---------------------!
    else

       cff  = one/real(nz, kind=rp)
       cff1 = one/sinh(theta_s)
       cff2 = hlf/tanh(hlf*theta_s)

       do i = 1-nh,nx+nh
          do j = 1-nh,ny+nh

             do k = 1,nz

                sc_w = cff*real(k-1-nz,kind=rp)
                sc_r = cff*(real(k-nz,kind=rp)-hlf)

                cs_r = (1.-theta_b)*cff1*sinh(theta_s*sc_r) &
                     +theta_b*(cff2*tanh(theta_s*(sc_r+hlf))-hlf)

                cs_w = (1.-theta_b)*cff1*sinh(theta_s*sc_w) &
                     +theta_b*(cff2*tanh(theta_s*(sc_w+hlf))-hlf)

                cff_w = hc * sc_w
                cff_r = hc * sc_r

                z_w0 = cff_w + cs_w*h(j,i) 
                z_r0 = cff_r + cs_r*h(j,i)  

                hinv = one / (h(j,i)+hc)

                zw(k,j,i) = z_w0 * (h(j,i)*hinv)
                zr(k,j,i) = z_r0 * (h(j,i)*hinv)

             enddo ! k

             zw(nz+1,j,i) = nul

          enddo ! j
       enddo ! i

    endif

  end subroutine setup_zr_zw_croco

end module mg_zr_zw
#else
        module mg_zr_zw_empty
        end module
#endif
