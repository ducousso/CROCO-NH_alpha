      subroutine rhs3d (tile)
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
      call rhs3d_tile (Istr,Iend,Jstr,Jend,
     &                                    A3d(1,1,trd), A3d(1,2,trd),
     &                      A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &                      A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &                      A2d(1,4,trd), A2d(1,5,trd), A2d(1,6,trd))
      return
      end
      subroutine rhs3d_tile (Istr,Iend,Jstr,Jend, ru,rv, CF,FC,DC,
     &                                 wrk1,wrk2, UFx,UFe, VFx,VFe)
      implicit none
      integer*4 Istr,Iend,Jstr,Jend, i,j,k
     &                   ,imin,imax,jmin,jmax
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
      real ru(Istr-2:Iend+2,Jstr-2:Jend+2,N),
     &     rv(Istr-2:Iend+2,Jstr-2:Jend+2,N),
     &     CF(Istr-2:Iend+2,0:N),  cff,cff1,cff2,
     &     FC(Istr-2:Iend+2,0:N),  cffX,
     &     DC(Istr-2:Iend+2,0:N),  cffE,
     &     wrk1(Istr-2:Iend+2,Jstr-2:Jend+2),    curvX,
     &     wrk2(Istr-2:Iend+2,Jstr-2:Jend+2),    curvE,
     &      UFx(Istr-2:Iend+2,Jstr-2:Jend+2),
     &      UFe(Istr-2:Iend+2,Jstr-2:Jend+2),
     &      VFx(Istr-2:Iend+2,Jstr-2:Jend+2),
     &      VFe(Istr-2:Iend+2,Jstr-2:Jend+2),    gamma
      parameter (gamma = -0.25D0)
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
      real rhoA(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rhoS(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /coup_rhoA/rhoA           /coup_rhoS/rhoS
      real rufrc(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rvfrc(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rufrc_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real rvfrc_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      common /coup_rufrc/rufrc         /coup_rvfrc/rvfrc
     &       /coup_rufrc_bak/rufrc_bak /coup_rvfrc_bak/rvfrc_bak
      real Zt_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real DU_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,5)
      real DV_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,5)
      real DU_avg2(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real DV_avg2(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /ocean_Zt_avg1/Zt_avg1
     &     /coup_DU_avg1/DU_avg1 /coup_DV_avg1/DV_avg1
     &     /coup_DU_avg2/DU_avg2 /coup_DV_avg2/DV_avg2
      real sustr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real svstr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real bustr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real bvstr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real bvstrg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen,    bms_tstart,   bms_tend,  tsbms,   
     &                                                            sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real stflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /forces_stflx/stflx
      real stflxg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2,NT)
      common /stfdat_stflxg/stflxg
      real stflxp(2,NT), stf_time(2,NT)
      real stf_cycle(NT), stf_scale(NT)
      integer*4 itstf(NT), stf_ncycle(NT), stf_rec(NT)
      integer*4 lstfgrd(NT), stf_tid(NT), stf_id(NT)
      common /stfdat1/ stflxp,  stf_time, stf_cycle, stf_scale
      common /stfdat2/ itstf, stf_ncycle, stf_rec, lstfgrd
      common /stfdat3/  stf_tid, stf_id
      real btflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /forces_btflx/btflx
      real srflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_srflx/srflx
      real srflxg(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      common /srfdat_srflxg/srflxg
      real srflxp(2),srf_time(2)
      real srf_cycle, srf_scale
      integer*4 itsrf, srf_ncycle, srf_rec
      integer*4 lsrfgrd, srf_tid, srf_id
      common /srfdat1/ srflxp, srf_time, srf_cycle, srf_scale
      common /srfdat2/ itsrf, srf_ncycle, srf_rec, lsrfgrd, srf_tid, 
     &                                                            srf_id
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
      do k=1,N
        if (WEST_INTER) then
          imin=IstrU-1
        else
          imin=max(IstrU-1,2)
        endif
        if (EAST_INTER) then
          imax=Iend+1
        else
          imax=min(Iend+1,Lmmpi)
        endif
        do j=Jstr,Jend
          do i=imin,imax
            wrk1(i,j)=(u(i-1,j,k,nrhs)-2.D0*u(i,j,k,nrhs)
     &               +u(i+1,j,k,nrhs))
            wrk2(i,j)=(Huon(i-1,j,k)-2.D0*Huon(i,j,k)
     &                +Huon(i+1,j,k))
          enddo
        enddo
        if (.not.WEST_INTER) then
          do j=Jstr,Jend
            wrk1(1,j)=wrk1(2,j)
            wrk2(1,j)=wrk2(2,j)
          enddo
        endif
        if (.not.EAST_INTER) then
          do j=Jstr,Jend
            wrk1(Lmmpi+1,j)=wrk1(Lmmpi,j)
            wrk2(Lmmpi+1,j)=wrk2(Lmmpi,j)
          enddo
        endif
        do j=Jstr,Jend
          do i=IstrU-1,Iend
            cffX=u(i,j,k,nrhs)+u(i+1,j,k,nrhs)
            if (cffX.gt.0.D0) then
              curvX=wrk1(i,j)
            else
              curvX=wrk1(i+1,j)
            endif
            UFx(i,j)=0.25D0*(cffX+gamma*curvX)*( Huon(i,j,k)+
     &            Huon(i+1,j,k)-0.125D0*(wrk2(i,j)+wrk2(i+1,j)))
          enddo
        enddo
        if (SOUTH_INTER) then
          jmin=JstrV-1
        else
          jmin=max(JstrV-1,2)
        endif
        if (NORTH_INTER) then
          jmax=Jend+1
        else
          jmax=min(Jend+1,Mmmpi)
        endif
        do j=jmin,jmax
          do i=Istr,Iend
            wrk1(i,j)=(v(i,j-1,k,nrhs)-2.D0*v(i,j,k,nrhs)
     &               +v(i,j+1,k,nrhs))
            wrk2(i,j)=(Hvom(i,j-1,k)-2.D0*Hvom(i,j,k)
     &                +Hvom(i,j+1,k))
          enddo
        enddo
        if (.not.SOUTH_INTER) then
          do i=Istr,Iend
            wrk1(i,1)=wrk1(i,2)
            wrk2(i,1)=wrk2(i,2)
          enddo
        endif
        if (.not.NORTH_INTER) then
          do i=Istr,Iend
            wrk1(i,Mmmpi+1)=wrk1(i,Mmmpi)
            wrk2(i,Mmmpi+1)=wrk2(i,Mmmpi)
          enddo
        endif
        do j=JstrV-1,Jend
          do i=Istr,Iend
            cffE=v(i,j,k,nrhs)+v(i,j+1,k,nrhs)
            if (cffE.gt.0.D0) then
              curvE=wrk1(i,j)
            else
              curvE=wrk1(i,j+1)
            endif
            VFe(i,j)=0.25D0*(cffE+gamma*curvE)*( Hvom(i,j,k)+
     &             Hvom(i,j+1,k)-0.125D0*(wrk2(i,j)+wrk2(i,j+1)))
          enddo
        enddo
        if (SOUTH_INTER) then
          jmin=Jstr-1
        else
          jmin=max(Jstr-1,1)
        endif
        if (NORTH_INTER) then
          jmax=Jend+1
        else
          jmax=min(Jend+1,Mmmpi)
        endif
        do j=jmin,jmax
          do i=IstrU,Iend
            wrk1(i,j)=(u(i,j-1,k,nrhs)-u(i,j,k,nrhs))
     &              +(u(i,j+1,k,nrhs)-u(i,j,k,nrhs))
          enddo
        enddo
        if (.not.SOUTH_INTER) then
          do i=IstrU,Iend
            wrk1(i,0)=wrk1(i,1)
          enddo
        endif
        if (.not.NORTH_INTER) then
          do i=IstrU,Iend
            wrk1(i,Mmmpi+1)=wrk1(i,Mmmpi)
          enddo
        endif
        do j=Jstr,Jend+1
          do i=IstrU-1,Iend
           wrk2(i,j)=Hvom(i-1,j,k)-2.D0*Hvom(i,j,k)+Hvom(i+1,j,k)
          enddo
        enddo
        do j=Jstr,Jend+1
          do i=IstrU,Iend
            cffX=u(i,j,k,nrhs)+u(i,j-1,k,nrhs)
            cffE=Hvom(i,j,k)+Hvom(i-1,j,k)
            if (cffE.gt.0.D0) then
              curvX=wrk1(i,j-1)
            else
              curvX=wrk1(i,j)
            endif
            UFe(i,j)=0.25D0*(cffX+gamma*curvX)*(cffE-0.125D0*(
     &                             wrk2(i,j)+wrk2(i-1,j) ))
          enddo
        enddo
        if (WEST_INTER) then
          imin=Istr-1
        else
          imin=max(Istr-1,1)
        endif
        if (EAST_INTER) then
          imax=Iend+1
        else
          imax=min(Iend+1,Lmmpi)
        endif
        do j=JstrV,Jend
          do i=imin,imax
            wrk1(i,j)=(v(i-1,j,k,nrhs)-v(i,j,k,nrhs))
     &              +(v(i+1,j,k,nrhs)-v(i,j,k,nrhs))
          enddo
        enddo
        if (.not.WEST_INTER) then
          do j=JstrV,Jend
            wrk1(0,j)=wrk1(1,j)
          enddo
        endif
        if (.not.EAST_INTER) then
          do j=JstrV,Jend
            wrk1(Lmmpi+1,j)=wrk1(Lmmpi,j)
          enddo
        endif
        do j=JstrV-1,Jend
          do i=Istr,Iend+1
           wrk2(i,j)=Huon(i,j-1,k)-2.D0*Huon(i,j,k)+Huon(i,j+1,k)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend+1
            cffE=v(i,j,k,nrhs)+v(i-1,j,k,nrhs)
            cffX=Huon(i,j,k)+Huon(i,j-1,k)
            if (cffX.gt.0.D0) then
              curvE=wrk1(i-1,j)
            else
              curvE=wrk1(i,j)
            endif
            VFx(i,j)=0.25D0*(cffE+gamma*curvE)*(cffX-0.125D0*(
     &                             wrk2(i,j)+wrk2(i,j-1) ))
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            ru(i,j,k)=ru(i,j,k)-UFx(i,j  )+UFx(i-1,j)
     &                         -UFe(i,j+1)+UFe(i  ,j)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            rv(i,j,k)=rv(i,j,k)-VFx(i+1,j)+VFx(i,j  )
     &                         -VFe(i  ,j)+VFe(i,j-1)
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do k=1,N
          do i=IstrU,Iend
            DC(i,k)=0.5625D0*(Hz(i  ,j,k)+Hz(i-1,j,k))
     &             -0.0625D0*(Hz(i+1,j,k)+Hz(i-2,j,k))
          enddo
        enddo
        do i=IstrU,Iend
          FC(i,0)=0.D0
          CF(i,0)=0.D0
        enddo
        do k=1,N-1,+1
          do i=IstrU,Iend
            cff=1.D0/(2.D0*DC(i,k+1)+DC(i,k)*(2.D0-FC(i,k-1)))
            FC(i,k)=cff*DC(i,k+1)
            CF(i,k)=cff*(6.D0*(u(i,j,k+1,nrhs)-u(i,j,k,nrhs))
     &                                  -DC(i,k)*CF(i,k-1))
          enddo
        enddo
        do i=IstrU,Iend
          CF(i,N)=0.D0
        enddo
        do k=N-1,1,-1
          do i=IstrU,Iend
            CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
          enddo
        enddo
        do k=1,N-1
          do i=IstrU,Iend
            FC(i,k)=( 0.5625D0*(We(i  ,j,k)+We(i-1,j,k))-0.0625D0*(
     &                                 We(i+1,j,k)+We(i-2,j,k) ))
     &             *( u(i,j,k,nrhs)+DC(i,k)*(
     &                              0.33333333333333D0*CF(i,k  )
     &                             +0.16666666666667D0*CF(i,k-1)
     &                                                       ))
          enddo
        enddo
        do i=IstrU,Iend
          FC(i,N)=0.D0
          FC(i,0)=0.D0
        enddo
        do k=1,N
          do i=IstrU,Iend
            ru(i,j,k)=ru(i,j,k)-FC(i,k)+FC(i,k-1)
          enddo
        enddo
        if (j.ge.JstrV) then
          do k=1,N
            do i=Istr,Iend
              DC(i,k)=0.5625D0*(Hz(i  ,j,k)+Hz(i,j-1,k))
     &               -0.0625D0*(Hz(i,j+1,k)+Hz(i,j-2,k))
            enddo
          enddo
          do i=Istr,Iend
            FC(i,0)=0.D0
            CF(i,0)=0.D0
          enddo
          do k=1,N-1,+1
            do i=Istr,Iend
              cff=1.D0/(2.D0*DC(i,k+1)+DC(i,k)*(2.D0-FC(i,k-1)))
              FC(i,k)=cff*DC(i,k+1)
              CF(i,k)=cff*(6.D0*(v(i,j,k+1,nrhs)-v(i,j,k,nrhs))
     &                                    -DC(i,k)*CF(i,k-1))
            enddo
          enddo
          do i=Istr,Iend
            CF(i,N)=0.D0
          enddo
          do k=N-1,1,-1
            do i=Istr,Iend
              CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
            enddo
          enddo
          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=( 0.5625D0*(We(i,j  ,k)+We(i,j-1,k))-0.0625D0*(
     &                                   We(i,j+1,k)+We(i,j-2,k) ))
     &               *( v(i,j,k,nrhs)+DC(i,k)*(
     &                                0.33333333333333D0*CF(i,k  )
     &                               +0.16666666666667D0*CF(i,k-1)
     &                                                         ))
            enddo
          enddo
          do i=Istr,Iend
            FC(i,N)=0.D0
            FC(i,0)=0.D0
          enddo
          do k=1,N
            do i=Istr,Iend
              rv(i,j,k)=rv(i,j,k)-FC(i,k)+FC(i,k-1)
            enddo
          enddo
        endif
        do i=IstrU,Iend
          rufrc(i,j)=ru(i,j,1)
     &        +(sustr(i,j)-bustr(i,j))*om_u(i,j)*on_u(i,j)
        enddo
        do k=2,N
          do i=IstrU,Iend
            rufrc(i,j)=rufrc(i,j)+ru(i,j,k)
          enddo
        enddo
        if (j.ge.JstrV) then
          do i=Istr,Iend
            rvfrc(i,j)=rv(i,j,1)
     &       +(svstr(i,j)-bvstr(i,j))*om_v(i,j)*on_v(i,j)
          enddo
          do k=2,N
            do i=Istr,Iend
              rvfrc(i,j)=rvfrc(i,j)+rv(i,j,k)
            enddo
          enddo
        endif
      enddo
      return
      end
