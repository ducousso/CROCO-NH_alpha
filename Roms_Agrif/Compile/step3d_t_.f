      subroutine step3d_t (tile)
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
      call step3d_t_tile (Istr,Iend,Jstr,Jend,
     &                    A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &                    A2d(1,4,trd), A2d(1,5,trd), A2d(1,6,trd),
     &                    A2d(1,7,trd), A2d(1,8,trd), A2d(1,9,trd),
     &                                                A3d(1,1,trd))
      return
      end
      subroutine step3d_t_tile (Istr,Iend,Jstr,Jend,
     &                          FX,FE, WORK, FC,CF,BC,DC,EC,GC, swdk)
      implicit none
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
      real visc2_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_sponge_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_sponge_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
      common /mixing_visc2_sponge_p/visc2_sponge_p
      real Akv(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Akt(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
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
      integer*4 Istr,Iend,Jstr,Jend, itrc, i,j,k, indx, kmld
     &       ,imin,imax,jmin,jmax,iAkt,nadv
      real FX(Istr-2:Iend+2,Jstr-2:Jend+2),
     &     FE(Istr-2:Iend+2,Jstr-2:Jend+2),   cff,
     &     WORK(Istr-2:Iend+2,Jstr-2:Jend+2), epsil,
     &     FC(Istr-2:Iend+2,0:N),
     &     CF(Istr-2:Iend+2,0:N),
     &     BC(Istr-2:Iend+2,0:N),
     &     DC(Istr-2:Iend+2,0:N),
     &     EC(Istr-2:Iend+2,0:N),
     &     GC(Istr-2:Iend+2,0:N),
     &     swdk(Istr-2:Iend+2,Jstr-2:Jend+2,0:N)
      real cff1,cff2,gama,dRz,hbltmp,sig,dXmax,dEmax,
     &     dpth,smax,amax,amaxx
      parameter (epsil=1.D-16)
      REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2
      REAL    :: ua, vel, cdiff, cdif
      REAL    :: flux1, flux2, flux3, flux4, flux5, flux6
      REAL    :: flx2, flx3, flx4, flx5
      REAL    :: mask0, mask1, mask2, mask3
      flux2(q_im1, q_i, ua, cdiff) = 0.5D0*( q_i + q_im1 )
      flux1(q_im1, q_i, ua, cdiff) = flux2(q_im1, q_i, ua, cdiff) -
     &      0.5D0*cdiff*sign(1.D0,ua)*(q_i-q_im1)
      flux4(q_im2, q_im1, q_i, q_ip1, ua) =
     &      ( 7.D0*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0D0
      flux3(q_im2, q_im1, q_i, q_ip1, ua) =
     &      flux4(q_im2, q_im1, q_i, q_ip1, ua) +
     &      sign(1.D0,ua)*((q_ip1 -
     &      q_im2)-3.D0*(q_i-q_im1))/12.0D0
      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      ( 37.D0*(q_i+q_im1) - 8.D0*(q_ip1+q_im2)
     &      +(q_ip2+q_im3) )/60.0D0
      flux5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
     &      -sign(1.D0,ua)*(
     &      (q_ip2-q_im3)-5.D0*(q_ip1-q_im2)+10.D0*(q_i-q_im1) )/60.0D0
      REAL    :: flux3_weno, flux5_weno
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
      nadv = 3
      do itrc=1,NT
        do k=1,N
            cdif=0.2D0
          if (SOUTH_INTER) then
            jmin=1
          else
            jmin=3
          endif
          if (NORTH_INTER) then
            jmax=Mmmpi+1
          else
            jmax=Mmmpi-1
          endif
          if (WEST_INTER) then
            imin=1
          else
            imin=3
          endif
          if (EAST_INTER) then
            imax=Lmmpi+1
          else
            imax=Lmmpi-1
          endif
          DO j = Jstr,Jend+1
            IF ( j.ge.jmin .and. j.le.jmax ) THEN
              DO i = Istr,Iend
                vel = Hvom(i,j,k)
                flx5 = vel*flux5_weno(
     &             t(i,j-3,k,nrhs,itrc), t(i,j-2,k,nrhs,itrc),
     &             t(i,j-1,k,nrhs,itrc), t(i,j  ,k,nrhs,itrc),
     &             t(i,j+1,k,nrhs,itrc), t(i,j+2,k,nrhs,itrc),  vel )
                FE(i,j)=flx5
              ENDDO
            ELSE IF ( j.eq.jmin-2 ) THEN
              DO i = Istr,Iend
                vel = Hvom(i,j,k)
                FE(i,j) = vel*flux1(
     &             t(i,j-1,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
              ENDDO
            ELSE IF ( j.eq.jmin-1 .and. jmax.ge.jmin ) THEN
              DO i = Istr,Iend
                vel = Hvom(i,j,k)
                flx3 = vel*flux3(
     &             t(i,j-2,k,nrhs,itrc), t(i,j-1,k,nrhs,itrc),
     &             t(i,j  ,k,nrhs,itrc), t(i,j+1,k,nrhs,itrc),  vel )
                FE(i,j)=flx3
              ENDDO
            ELSE IF ( j.eq.jmax+2 ) THEN
              DO i = Istr,Iend
                vel = Hvom(i,j,k)
                FE(i,j) = vel*flux1(
     &             t(i,j-1,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
              ENDDO
            ELSE IF ( j.eq.jmax+1 ) THEN
              DO i = Istr,Iend
                vel = Hvom(i,j,k)
                flx3 = vel*flux3(
     &             t(i,j-2,k,nrhs,itrc), t(i,j-1,k,nrhs,itrc),
     &             t(i,j  ,k,nrhs,itrc), t(i,j+1,k,nrhs,itrc),  vel )
                FE(i,j)=flx3
              ENDDO
            ENDIF
          ENDDO
          DO i = Istr,Iend+1
            IF ( i.ge.imin .and. i.le.imax ) THEN
              DO j = Jstr,Jend
                vel = Huon(i,j,k)
                flx5 = vel*flux5_weno(
     &             t(i-3,j,k,nrhs,itrc), t(i-2,j,k,nrhs,itrc),
     &             t(i-1,j,k,nrhs,itrc), t(i  ,j,k,nrhs,itrc),
     &             t(i+1,j,k,nrhs,itrc), t(i+2,j,k,nrhs,itrc),  vel )
                FX(i,j)=flx5
              ENDDO
            ELSE IF ( i.eq.imin-2 ) THEN
              DO j = Jstr,Jend
                vel = Huon(i,j,k)
                FX(i,j) = vel*flux1(
     &             t(i-1,j,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
              ENDDO
            ELSE IF ( i.eq.imin-1 .and. imax.ge.imin ) THEN
              DO j = Jstr,Jend
                vel = Huon(i,j,k)
                flx3 = vel*flux3(
     &             t(i-2,j,k,nrhs,itrc), t(i-1,j,k,nrhs,itrc),
     &             t(i  ,j,k,nrhs,itrc), t(i+1,j,k,nrhs,itrc),  vel )
                FX(i,j)=flx3
              ENDDO
            ELSE IF ( i.eq.imax+2 ) THEN
              DO j = Jstr,Jend
                vel = Huon(i,j,k)
                FX(i,j) = vel*flux1(
     &             t(i-1,j,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
              ENDDO
            ELSE IF ( i.eq.imax+1 ) THEN
              DO j = Jstr,Jend
                vel = Huon(i,j,k)
                flx3 = vel*flux3(
     &             t(i-2,j,k,nrhs,itrc), t(i-1,j,k,nrhs,itrc),
     &             t(i  ,j,k,nrhs,itrc), t(i+1,j,k,nrhs,itrc),  vel )
                FX(i,j)=flx3
              ENDDO
            ENDIF
          ENDDO
          do j=Jstr,Jend
            do i=Istr,Iend
              t(i,j,k,nnew,itrc)=Hz_bak(i,j,k)*t(i,j,k,nstp,itrc)
     &                     -dt*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                          +FE(i,j+1)-FE(i,j)
     &                                                           )
            enddo
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do itrc=1,NT
          do k=3,N-3
            do i=Istr,Iend
              FC(i,k)=We(i,j,k)*
     &              flux5_weno(
     &             t(i,j,k-2,nadv,itrc), t(i,j,k-1,nadv,itrc),
     &             t(i,j,k  ,nadv,itrc), t(i,j,k+1,nadv,itrc),
     &             t(i,j,k+2,nadv,itrc), t(i,j,k+3,nadv,itrc), 
     &                                                        We(i,j,k))
            enddo
          enddo
          do i=Istr,Iend
            FC(i,2)=We(i,j,2)*
     &               flux3_weno(
     &           t(i,j,1,nadv,itrc), t(i,j,2,nadv,itrc),
     &           t(i,j,3,nadv,itrc), t(i,j,4,nadv,itrc), We(i,j,2))
            FC(i,N-2)=We(i,j,N-2)*
     &               flux3_weno(
     &           t(i,j,N-3,nadv,itrc), t(i,j,N-2,nadv,itrc),
     &           t(i,j,N-1,nadv,itrc), t(i,j,N  ,nadv,itrc), 
     &                                                      We(i,j,N-2))
            FC(i,  1)=0.5D0*We(i,j,  1)*( t(i,j,1  ,nadv,itrc)
     &                               +  t(i,j,2  ,nadv,itrc))
            FC(i,N-1)=0.5D0*We(i,j,N-1)*( t(i,j,N-1,nadv,itrc)
     &                               +  t(i,j,N,  nadv,itrc))
            FC(i,0)=0.D0
            FC(i,N )=0.D0
            CF(i,0)=dt*pm(i,j)*pn(i,j)
          enddo
          do k=1,N
            do i=Istr,Iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-CF(i,0)*(FC(i,k)
     &                                                      -FC(i,k-1))
            enddo
          enddo
          do i=Istr,Iend
            FC(i,N)=dt*stflx(i,j,itrc)
            FC(i,0)=-dt*btflx(i,j,itrc)
          enddo
          if (itrc.eq.itemp) then
            do k=1,N-1
              do i=Istr,Iend
                FC(i,k)=0.D0
              enddo
            enddo
          endif
          if (itrc.eq.itemp) then
            do k=1,N
              do i=Istr,Iend
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+FC(i,k )
     &                                               -FC(i,k-1)
              enddo
            enddo
          endif
          indx=min(itrc,itemp)
          do i=istr,iend
            FC(i,1)=dt* Akt(i,j,1,indx)
     &                               /( z_r(i,j,2)-z_r(i,j,1) )
            cff=1.D0/(Hz(i,j,1)+FC(i,1))
            CF(i,1)= cff*FC(i,1)
            DC(i,1)= cff* t(i,j,1,nnew,itrc)
          enddo
          do k=2,N-1,+1
            do i=istr,iend
              FC(i,k)=dt* Akt(i,j,k,indx)
     &                              /( z_r(i,j,k+1)-z_r(i,j,k) )
              cff=1.D0/(Hz(i,j,k)+FC(i,k)+FC(i,k-1)*(1.D0-CF(i,k-1)))
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(t(i,j,k,nnew,itrc)+FC(i,k-1)*DC(i,k-1)
     &                                                          )
            enddo
          enddo
          do i=istr,iend
             t(i,j,N,nnew,itrc)=( t(i,j,N,nnew,itrc)
     &                           +FC(i,N-1)*DC(i,N-1) )
     &                        /(Hz(i,j,N)+FC(i,N-1)*(1.D0-CF(i,N-1)))
        enddo
          do k=N-1,1,-1
            do i=istr,iend
              t(i,j,k,nnew,itrc)=DC(i,k)+CF(i,k)*t(i,j,k+1,nnew,itrc)
            enddo
          enddo
        enddo
      enddo
      do itrc=1,NT
        call t3dbc_tile (Istr,Iend,Jstr,Jend, nnew,itrc, WORK)
      enddo
      do itrc=1,NT
        call exchange_r3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                          t(-2,-2,1,nnew,itrc))
      enddo
      return
      end
      function flux5_weno(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
      implicit none
      REAL    :: flux5_weno
      REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua
      REAL    :: IS0, IS1, IS2
      REAl    :: d0, d1, d2
      REAl    :: a0, a1, a2
      REAL    :: w0, w1, w2
      REAL    :: p0, p1, p2
      REAL    :: Eps, cff1, cff2, T5
      Eps = 1.D-40
      d0=1.D0/10.D0
      d1=6.D0/10.D0
      d2=3.D0/10.D0
      cff1=13.D0/12.D0
      cff2=1.D0/6.D0
      if (ua .ge. 0.D0) then
        IS0 = cff1*(q_im3 - 2.D0*q_im2 + q_im1)**2
     &                        + 0.25D0*(q_im3 - 4.D0*q_im2 + 3*q_im1)**2
        IS1 = cff1*(q_im2 - 2.D0*q_im1 + q_i)**2
     &                        + 0.25D0*(q_im2 - q_i)**2
        IS2 = cff1*(q_im1 - 2.D0*q_i + q_ip1)**2
     &                        + 0.25D0*(3.D0*q_im1 - 4.D0*q_i + 
     &                                                         q_ip1)**2
        T5 = abs(IS2-IS0)
        a0  = d0*(1+T5/(Eps+IS0))
        a1  = d1*(1+T5/(Eps+IS1))
        a2  = d2*(1+T5/(Eps+IS2))
        w0  = a0/(a0+a1+a2)
        w1  = a1/(a0+a1+a2)
        w2  = a2/(a0+a1+a2)
        p0  = cff2*(2.D0*q_im3 - 7.D0*q_im2 + 11.D0*q_im1)
        p1  = cff2*(-q_im2   + 5.D0*q_im1 + 2.D0*q_i)
        p2  = cff2*(2.D0*q_im1 + 5.D0*q_i   - q_ip1)
        flux5_weno = w0*p0 + w1*p1 +w2*p2
      else
        IS0 = cff1*(q_ip2 - 2.D0*q_ip1 + q_i)**2
     &                         + 0.25D0*(q_ip2 -4.D0*q_ip1 + 3*q_i)**2
        IS1 = cff1*(q_ip1 - 2.D0*q_i + q_im1)**2
     &                         + 0.25D0*(q_ip1-q_im1)**2
        IS2 = cff1*(q_i   - 2.D0*q_im1 + q_im2)**2
     &                         + 0.25D0*(3.D0*q_i -4.D0*q_im1 + 
     &                                                         q_im2)**2
        T5 = abs(IS2-IS0)
        a0  = d0*(1+T5/(Eps+IS0))
        a1  = d1*(1+T5/(Eps+IS1))
        a2  = d2*(1+T5/(Eps+IS2))
        w0  = a0/(a0+a1+a2)
        w1  = a1/(a0+a1+a2)
        w2  = a2/(a0+a1+a2)
        p0  = cff2*(2.D0*q_ip2 - 7.D0*q_ip1 + 11.D0*q_i)
        p1  = cff2*(-q_ip1   + 5.D0*q_i   + 2.D0*q_im1)
        p2  = cff2*(2.D0*q_i   + 5.D0*q_im1 - q_im2)
        flux5_weno = w0*p0 + w1*p1 +w2*p2
      endif
      return
      end
      function flux3_weno( q_im2, q_im1, q_i, q_ip1, ua)
      implicit none
      REAL    :: flux3_weno
      REAL    :: q_im2, q_im1, q_i, q_ip1, ua
      REAL    :: IS0, IS1
      REAl    :: a0, a1
      REAL    :: w0, w1
      REAL    :: p0, p1
      REAL    :: Eps, d0,d1, T3
      Eps = 1.D-40
      d0=1.D0/3.D0
      d1=2.D0/3.D0
      if (ua .ge. 0.D0) then
        IS0 = (q_im1-q_im2)**2
        IS1 = (q_im1-q_i)**2
        T3 = abs(IS1-IS0)
        a0  = d0*(1+T3/(Eps+IS0))
        a1  = d1*(1+T3/(Eps+IS1))
        w0  = a0/(a0+a1)
        w1  = a1/(a0+a1)
        p0  = 1.D0/2.D0*(3.D0*q_im1-q_im2)
        p1  = 1.D0/2.D0*(q_im1+q_i)
        flux3_weno = w0*p0 + w1*p1
      else
        IS0 = (q_i-q_ip1)**2
        IS1 = (q_im1-q_i)**2
        T3 = abs(IS1-IS0)
        a0  = d0*(1+T3/(Eps+IS0))
        a1  = d1*(1+T3/(Eps+IS1))
        w0  = a0/(a0+a1)
        w1  = a1/(a0+a1)
        p0  = 1.D0/2.D0*(3.D0*q_i-q_ip1)
        p1  = 1.D0/2.D0*(q_im1+q_i)
        flux3_weno = w0*p0 + w1*p1
      endif
      return
      end
