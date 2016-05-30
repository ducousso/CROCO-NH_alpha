      subroutine prsgrd (tile)
      implicit none
      integer*4 tile, trd, omp_get_thread_num
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=128,  MMm0=4,    N=40)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=1, NP_ETA=1,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6, Np=N+1)
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, itemp
      integer*4   ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      parameter (itemp=1)
      parameter (ntrc_salt=0)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_sed=0)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,5,0:NPP-1)
      common /private_scratch/ A2d,A3d
      integer*4 chunk_size_X,margin_X,chunk_size_E,margin_E
      integer*4 Istr,Iend,Jstr,Jend, i_X,j_E
      chunk_size_X=(Lmmpi+NSUB_X-1)/NSUB_X
      margin_X=(NSUB_X*chunk_size_X-Lmmpi)/2
      chunk_size_E=(Mmmpi+NSUB_E-1)/NSUB_E
      margin_E=(NSUB_E*chunk_size_E-Mmmpi)/2
      j_E=tile/NSUB_X
      i_X=tile-j_E*NSUB_X
      Istr=1+i_X*chunk_size_X-margin_X
      Iend=Istr+chunk_size_X-1
      Istr=max(Istr,1)
      Iend=min(Iend,Lmmpi)
      Jstr=1+j_E*chunk_size_E-margin_E
      Jend=Jstr+chunk_size_E-1
      Jstr=max(Jstr,1)
      Jend=min(Jend,Mmmpi)
      trd=omp_get_thread_num()
      call prsgrd_tile (Istr,Iend,Jstr,Jend,
     &                  A3d(1,1,trd), A3d(1,2,trd),
     &                  A2d(1,1,trd), A2d(1,2,trd)
     &                                           )
      return
      end
      subroutine prsgrd_tile (Istr,Iend,Jstr,Jend,
     &                          ru,rv,rsurf, pgrd
     &                                           )
      implicit none
      integer*4 Istr,Iend,Jstr,Jend, i,j,k
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=128,  MMm0=4,    N=40)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=1, NP_ETA=1,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6, Np=N+1)
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, itemp
      integer*4   ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      parameter (itemp=1)
      parameter (ntrc_salt=0)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_sed=0)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      real ru(Istr-2:Iend+2,Jstr-2:Jend+2,N),  cff,
     &     rv(Istr-2:Iend+2,Jstr-2:Jend+2,N),  cff1,
     &  rsurf(Istr-2:Iend+2,Jstr-2:Jend+2),    gamma,
     &   pgrd(Istr-2:Iend+2)
      real h(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real hinv(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real f(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real fomn(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real xp(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real xr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real yp(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real yr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /grid_xr/xr /grid_xp/xp /grid_yp/yp /grid_yr/yr
      real pm(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pm_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pm_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v
      real pmon_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmon_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmon_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real grdscl(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl
      real u(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real v(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real t(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hz_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_w(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Huon(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hvom(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real We(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real rho1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real rho(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real dt, dtfast, time, time2, time_start, tdays
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &                       ndtfast, iic, kstp, krhs, knew, next_kstp,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zob
     &               , xl, el, visc2, visc4, gamma2
      real  theta_s,   theta_b,   Tcline,  hc
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
      real  rx0, rx1
      real  tnu2(NT),tnu4(NT)
      real R0,T0,S0, Tcoef, Scoef
      real weight(6,0:NWEIGHT)
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
      logical ldefhis
      logical got_tini(NT)
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zob,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &           , theta_s,   theta_b,   Tcline,  hc
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1,       tnu2,    tnu4
     &                      , R0,T0,S0,  Tcoef,   Scoef
     &                      , weight
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                      , got_tini
     &                      , ldefhis
      real Akv_bak
      real Akt_bak(NT)
      common /scalars_akt/ Akv_bak, Akt_bak
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real*8 volume, avgke, avgpe, avgkp, bc_crss
     &        , avg_vol, avg_rho
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
     &        , avg_vol, avg_rho
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers_roms/CPU_time,proc,trd_count
      logical EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      integer*4 mynode, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE, p_NW,p_NE
      common /comm_setup/ mynode, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE,
     &  p_NW,p_NE, EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  day2sec=86400.D0,
     &           sec2day=1.D0/86400.D0, jul_off=2440000.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-9999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      integer*4 IstrR,IendR,JstrR,JendR
      integer*4 IstrU
      integer*4 JstrV
      if (.not.WEST_INTER) then
        IstrR=Istr-1
        IstrU=Istr+1
      else
        IstrR=Istr
        IstrU=Istr
      endif
      if (.not.EAST_INTER) then
        IendR=Iend+1
      else
        IendR=Iend
      endif
      if (.not.SOUTH_INTER) then
        JstrR=Jstr-1
        JstrV=Jstr+1
      else
        JstrR=Jstr
        JstrV=Jstr
      endif
      if (.not.NORTH_INTER) then
        JendR=Jend+1
      else
        JendR=Jend
      endif
      do j=JstrV-1,Jend
        do i=IstrU-1,Iend
          rsurf(i,j)=rho(i,j,N) + (rho(i,j,N)-rho(i,j,N-1))
     &                             *(z_w(i,j,N)-z_r(i,j,N))
     &                           /(z_r(i,j,N)-z_r(i,j,N-1))
        enddo
        if (j.ge.Jstr) then
          cff=0.5D0*g/rho0
          do i=IstrU,Iend
            pgrd(i)=(g+cff*(rsurf(i-1,j)+rsurf(i,j)))*( z_w(i-1,j,N)
     &                                                   -z_w(i,j,N))
     &     +cff*( (rho(i-1,j,N)-rsurf(i,j))*(z_w(i-1,j,N)-z_r(i,j,N))
     &           +(rsurf(i-1,j)-rho(i,j,N))*(z_w(i,j,N)-z_r(i-1,j,N))
     &                                                              )
            ru(i,j,N)=0.5D0*(Hz(i,j,N)+Hz(i-1,j,N))*on_u(i,j)*pgrd(i)
          enddo
          cff=0.5D0*g/rho0
          do k=N-1,1,-1
            do i=IstrU,Iend
              pgrd(i)=pgrd(i)-cff*(
     &                                   (rho(i,j,k+1)-rho(i-1,j,k))
     &                                  *(z_r(i-1,j,k+1)-z_r(i,j,k))
     &                                  +(rho(i,j,k)-rho(i-1,j,k+1))
     &                                  *(z_r(i,j,k+1)-z_r(i-1,j,k))
     &                                                             )
              ru(i,j,k)=0.5D0*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)*pgrd(i)
            enddo
          enddo
        endif
        if (j.ge.JstrV) then
          cff=0.5D0*g/rho0
          do i=Istr,Iend
            pgrd(i)=(g+cff*(rsurf(i,j-1)+rsurf(i,j)))*( z_w(i,j-1,N)
     &                                                   -z_w(i,j,N))
     &     +cff*( (rho(i,j-1,N)-rsurf(i,j))*(z_w(i,j-1,N)-z_r(i,j,N))
     &           +(rsurf(i,j-1)-rho(i,j,N))*(z_w(i,j,N)-z_r(i,j-1,N))
     &                                                              )
            rv(i,j,N)=0.5D0*(Hz(i,j,N)+Hz(i,j-1,N))*om_v(i,j)*pgrd(i)
          enddo
          cff=0.5D0*g/rho0
          do k=N-1,1,-1
            do i=Istr,Iend
              pgrd(i)=pgrd(i)-cff*(
     &                                   (rho(i,j,k+1)-rho(i,j-1,k))
     &                                  *(z_r(i,j-1,k+1)-z_r(i,j,k))
     &                                  +(rho(i,j,k)-rho(i,j-1,k+1))
     &                                  *(z_r(i,j,k+1)-z_r(i,j-1,k))
     &                                                             )
              rv(i,j,k)=0.5D0*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)*pgrd(i)
            enddo
          enddo
        endif
      enddo
      return
      end
