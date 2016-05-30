      subroutine init_scalars (ierr)
      implicit none
      integer*4 ierr, i,j,itrc, lvar, lenstr
      character*20 nametrc, unitt
      character*60 vname1, vname3
      integer*4 omp_get_num_threads
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
      integer*4 filetype_his, filetype_avg
     &       ,filetype_dia, filetype_dia_avg
     &       ,filetype_diaM, filetype_diaM_avg
     &       ,filetype_diabio, filetype_diabio_avg
      parameter (filetype_his=1, filetype_avg=2,
     &           filetype_dia=3, filetype_dia_avg=4,
     &           filetype_diaM=5, filetype_diaM_avg=6,
     &           filetype_diabio=7,filetype_diabio_avg=8)
      integer*4 iloop, indextemp
      integer*4 indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
      integer*4 indxU, indxV, indxT
      parameter (indxU=6, indxV=7, indxT=8)
      integer*4 indxBSD, indxBSS
      parameter (indxBSD=indxT+ntrc_salt+ntrc_pas+ntrc_bio+1,
     &           indxBSS=101)
      integer*4 indxO, indxW, indxR, indxVisc, indxDiff, indxAkv, 
     &                                                           indxAkt
      parameter (indxO=indxT+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed
     &                      +ntrc_diats+ntrc_diauv+ntrc_diabio+1,
     &           indxW=indxO+1, indxR=indxO+2, indxVisc=indxO+3,
     &           indxDiff=indxO+4,indxAkv=indxO+5, indxAkt=indxO+6)
      integer*4 indxSSH
      parameter (indxSSH=indxAkt+12)
      integer*4 indxSUSTR, indxSVSTR
      parameter (indxSUSTR=indxSSH+1, indxSVSTR=indxSSH+2)
      integer*4 indxTime2
      parameter (indxTime2=indxSSH+3)
      integer*4 indxShflx, indxShflx_rsw
      parameter (indxShflx=indxSSH+4)
      parameter (indxShflx_rsw=indxShflx+1)
      integer*4 indxSST, indxdQdSST
      parameter (indxSST=indxShflx_rsw+1, indxdQdSST=indxShflx_rsw+2)
      integer*4 indxWstr
      parameter (indxWstr=indxSUSTR+21)
      integer*4 indxUWstr
      parameter (indxUWstr=indxSUSTR+22)
      integer*4 indxVWstr
      parameter (indxVWstr=indxSUSTR+23)
      integer*4 indxBostr
      parameter (indxBostr=indxSUSTR+24)
      integer*4 indxWWA,indxWWD,indxWWP,indxWEB
      parameter (indxWWA=indxSUSTR+32, indxWWD=indxWWA+1,
     &           indxWWP=indxWWA+2
     &                             )
      integer*4 r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar, b3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,b3dvar=12)
      integer*4 xi_rho,xi_u, eta_rho,eta_v
      parameter (xi_rho=LLm+2,  xi_u=xi_rho-1,
     &           eta_rho=MMm+2, eta_v=eta_rho-1)
      integer*4 ncidfrc, ncidbulk, ncidclm,  ntsms ,
     &        ntsrf,  ntssh,  ntsst, ntsss, ntuclm,
     &        ntbulk, ncidqbar, ntqbar, ntww
      integer*4 nttclm(NT), ntstf(NT), nttsrc(NT)
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
      integer*4 rstT(NT)
      integer*4  ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisU,   hisV,   hisR,    hisHbl, hisHbbl
     &      , hisO,   hisW,   hisVisc, hisDiff
     &      , hisAkv, hisAkt, hisAks
      integer*4 hisT(NT)
      logical wrthis(500+NT)
      common/incscrum/
     &        ncidfrc, ncidbulk,ncidclm, ntsms, ntsrf, ntssh, ntsst
     &      , ntuclm, ntsss, ntbulk, ncidqbar, ntqbar, ntww
     &                        ,  nttclm, ntstf, nttsrc
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV,   rstT
     &      , ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisU,    hisV,     hisT,    hisR
     &      , hisO,    hisW,     hisVisc, hisDiff
     &      , hisAkv,  hisAkt,   hisAks
     &      , hisHbl,  hisHbbl
     &      , wrthis
      character*80 date_str, title, start_date
      character*80 ininame,  grdname,  hisname
     &         ,   rstname,  frcname,  bulkname,  usrname
     &         ,   qbarname, tsrcname
      character*75  vname(20, 500)
      common /cncscrum/       date_str,   title,  start_date
     &         ,   ininame,  grdname, hisname
     &         ,   rstname,  frcname, bulkname,  usrname
     &         ,   qbarname, tsrcname
     &                      ,  vname
      real bry_time(2)
      common /bry_indices_array/ bry_time
      real bry_cycle
      common /bry_indices_real/ bry_cycle
      integer*4 bry_id, bry_time_id, bry_ncycle, bry_rec, itbry, ntbry
      common /bry_indices_integer/ bry_id, bry_time_id, bry_ncycle,
     &                             bry_rec, itbry, ntbry
C$OMP PARALLEL
C$OMP CRITICAL (isca_cr_rgn)
      numthreads=omp_get_num_threads()
C$OMP END CRITICAL (isca_cr_rgn)
C$OMP END PARALLEL
      if (mynode.eq.0) write(stdout,'(1x,A,3(1x,A,I3),A)') 'NUMBER',
     &    'OF THREADS:',numthreads,'BLOCKING:',NSUB_X,'x',NSUB_E,'.'
      if (numthreads.gt.NPP) then
        if (mynode.eq.0) write(stdout,'(/1x,A,I3/)')
     &    'ERROR: Requested number of threads exceeds NPP =', NPP
        ierr=ierr+1
      elseif (mod(NSUB_X*NSUB_E,numthreads).ne.0) then
        if (mynode.eq.0) write(stdout,
     &                '(/1x,A,1x,A,I3,4x,A,I3,4x,A,I4,A)') 'ERROR:',
     &                'wrong choice of numthreads =', numthreads,
     &                'NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E, '.'
        ierr=ierr+1
      endif
      time=0.D0
      tdays=0.D0
      PREDICTOR_2D_STEP=.FALSE.
      iic=0
      kstp=1
      krhs=1
      knew=1
      ntstart=1
      nstp=1
      nrhs=1
      nnew=1
      nfast=1
      synchro_flag=.true.
      first_time=0
      may_day_flag=0
      do j=0,NPP
        do i=0,31
          proc(i,j)=0
          CPU_time(i,j)=0.D0
        enddo
      enddo
      trd_count=0
      nrecrst=0
      nrechis=0
      tile_count=0
      bc_count=0
      avgke=0.D0
      avgpe=0.D0
      avgkp=0.D0
      volume=0.D0
      hmin=+1.D+20
      hmax=-1.D+20
      grdmin=+1.D+20
      grdmax=-1.D+20
      Cu_min=+1.D+20
      Cu_max=-1.D+20
      bc_crss=0.D0
      rx0=-1.D+20
      rx1=-1.D+20
      ncidrst=-1
      ncidhis=-1
      ncidfrc=-1
      ncidbulk=-1
      ncidclm=-1
      ncidqbar=-1
       bry_id=-1
      call get_date (date_str)
      vname(:,:)='  '
      vname(1,indxTime)='scrum_time'
      vname(2,indxTime)='time since initialization'
      vname(3,indxTime)='second'
      vname(4,indxTime)='time, scalar, series'
      vname(5,indxTime)='time'
      vname(6,indxTime)='    '
      vname(7,indxTime)='T'
      vname(8,indxTime)='.NO.'
      vname(9,indxTime)=' '
      vname(1,indxTime2)='time'
      vname(2,indxTime2)='time since initialization'
      vname(3,indxTime2)='second'
      vname(4,indxTime2)='time, scalar, series'
      vname(5,indxTime2)='time'
      vname(6,indxTime2)='  '
      vname(7,indxTime2)='T '
      vname(8,indxTime2)='.NO.'
      vname(9,indxTime2)=' '
      vname(1,indxZ)='zeta                                        '
      vname(2,indxZ)='free-surface                                '
      vname(3,indxZ)='meter                                       '
      vname(4,indxZ)='free-surface, scalar, series                '
      vname(5,indxZ)='sea_surface_height                          '
      vname(6,indxZ)='x_rho y_rho                                 '
      vname(7,indxZ)='                                            '
      vname(1,indxUb)='ubar                                       '
      vname(2,indxUb)='vertically integrated u-momentum component '
      vname(3,indxUb)='meter second-1                             '
      vname(4,indxUb)='ubar-velocity, scalar, series              '
      vname(5,indxUb)='barotropic_sea_water_x_'
     &               // 'velocity_at_u_location'
      vname(7,indxUb)='                                           '
      vname(1,indxVb)='vbar                                       '
      vname(2,indxVb)='vertically integrated v-momentum component '
      vname(3,indxVb)='meter second-1                             '
      vname(4,indxVb)='vbar-velocity, scalar, series              '
      vname(5,indxVb)='barotropic_sea_water_y_velocity_at_'
     &                // 'v_location'
      vname(7,indxVb)='                                           '
      vname(11,indxVb)=' '
      vname(1,indxBostr)='bostr                                   '
      vname(2,indxBostr)='Kinematic bottom stress                 '
      vname(3,indxBostr)='N/m2                                    '
      vname(4,indxBostr)='                                        '
      vname(5,indxBostr)='                                        '
      vname(6,indxBostr)='x_rho y_rho                             '
      vname(7,indxBostr)='                                        '
      vname(1,indxWstr)='wstr                                     '
      vname(2,indxWstr)='Kinematic wind stress                    '
      vname(3,indxWstr)='N/m2                                     '
      vname(4,indxWstr)='                                         '
      vname(5,indxWstr)='magnitude_of_surface_downward_stress     '
      vname(6,indxWstr)='x_rho y_rho                              '
      vname(7,indxWstr)='                                         '
      vname(1,indxUWstr)='sustr                                   '
      vname(2,indxUWstr)='Kinematic u wind stress component       '
      vname(3,indxUWstr)='N/m2                                    '
      vname(4,indxUWstr)='                                        '
      vname(5,indxUWstr)='surface_downward_eastward_stress        '
      vname(7,indxUWstr)='                                        '
      vname(8,indxUWstr)='                                        '
      vname(1,indxVWstr)='svstr                                   '
      vname(2,indxVWstr)='Kinematic v wind stress component       '
      vname(3,indxVWstr)='N/m2                                    '
      vname(4,indxVWstr)='                                        '
      vname(5,indxVWstr)='surface_downward_northward_stress       '
      vname(7,indxVWstr)='                                        '
      vname(1,indxU)='u                                           '
      vname(2,indxU)='u-momentum component                        '
      vname(3,indxU)='meter second-1                              '
      vname(4,indxU)='u-velocity, scalar, series                  '
      vname(5,indxU)='sea_water_x_velocity_at_u_location          '
      vname(7,indxU)='                                            '
      vname(1,indxV)='v                                           '
      vname(2,indxV)='v-momentum component                        '
      vname(3,indxV)='meter second-1                              '
      vname(4,indxV)='v-velocity, scalar, series                  '
      vname(5,indxV)='sea_water_y_velocity_at_v_location          '
      vname(7,indxV)='                                            '
      vname(1,indxT)='temp                                        '
      vname(2,indxT)='potential temperature                       '
      vname(3,indxT)='Celsius                                     '
      vname(4,indxT)='temperature, scalar, series                 '
      vname(5,indxT)='sea_water_potential_temperature             '
      vname(6,indxT)='x_rho y_rho                                 '
      vname(7,indxT)='                                            '
      vname(1,indxShflx)='shflux                                  '
      vname(2,indxShflx)='surface net heat flux                   '
      vname(3,indxShflx)='Watts meter-2                           '
      vname(4,indxShflx)='surface heat flux, scalar, series       '
      vname(6,indxShflx)='x_rho y_rho                             '
      vname(7,indxShflx)='                                        '
      vname(1,indxShflx_rsw)='swrad                               '
      vname(2,indxShflx_rsw)='Short-wave surface radiation        '
      vname(3,indxShflx_rsw)='Watts meter-2                       '
      vname(4,indxShflx_rsw)='                                    '
      vname(5,indxShflx_rsw)='                                    '
      vname(6,indxShflx_rsw)='lat_rho lon_rho                     '
      vname(7,indxShflx_rsw)='                                    '
      vname(1,indxO)='omega                                       '
      vname(2,indxO)='S-coordinate vertical momentum component    '
      vname(3,indxO)='meter second-1                              '
      vname(4,indxO)='omega, scalar, series                       '
      vname(5,indxO)='                                            '
      vname(6,indxO)='lat_rho lon_rho                             '
      vname(7,indxO)='                                            '
      vname(1,indxW)='w                                           '
      vname(2,indxW)='vertical momentum component                 '
      vname(3,indxW)='meter second-1                              '
      vname(4,indxW)='w-velocity, scalar, series                  '
      vname(5,indxW)='upward_sea_water_velocity                   '
      vname(6,indxW)='lat_rho lon_rho                             '
      vname(1,indxR)='rho                                         '
      vname(2,indxR)='density anomaly                             '
      vname(3,indxR)='kilogram meter-3                            '
      vname(4,indxR)='density, scalar, series                     '
      vname(5,indxR)='sea_water_sigma_t                           '
      vname(7,indxR)='                                            '
      vname(1,indxAkv)='AKv                                       '
      vname(2,indxAkv)='vertical viscosity coefficient            '
      vname(3,indxAkv)='meter2 second-1                           '
      vname(4,indxAkv)='AKv, scalar, series                       '
      vname(5,indxAkv)='ocean_vertical_momentum_diffusivity_'
     &                 / / 'at_w_location                         '
      vname(6,indxAkv)='lat_rho lon_rho                           '
      vname(7,indxAkv)='                                          '
      vname(1,indxAkt)='AKt                                       '
      vname(2,indxAkt)='temperature vertical diffusion coefficient'
      vname(3,indxAkt)='meter2 second-1                           '
      vname(4,indxAkt)='AKt, scalar, series                       '
      vname(5,indxAkt)='ocean_vertical_heat_diffusivity_'
     &                               / /  'at_w_location          '
      vname(6,indxAkt)='lat_rho lon_rho                           '
      vname(7,indxAkt)='                                          '
      vname(1,indxSSH)='SSH                                       '
      vname(2,indxSSH)='sea surface height                        '
      vname(3,indxSSH)='meter                                     '
      vname(4,indxSSH)='SSH, scalar, series                       '
      vname(5,indxSSH)='sea_surface_height_above_sea_level        '
      vname(6,indxSSH)='lat_rho lon_rho                           '
      vname(7,indxSSH)='                                          '
      vname(1,indxSUSTR)='sustr                                   '
      vname(2,indxSUSTR)='surface u-momentum stress               '
      vname(3,indxSUSTR)='Newton meter-2                          '
      vname(4,indxSUSTR)='surface u-mom. stress, scalar, series   '
      vname(5,indxSUSTR)='surface_downward_x_stress               '
      vname(6,indxSUSTR)='lat_u lon_u                             '
      vname(7,indxSUSTR)='                                        '
      vname(1,indxSVSTR)='svstr                                   '
      vname(2,indxSVSTR)='surface v-momentum stress               '
      vname(3,indxSVSTR)='Newton meter-2                          '
      vname(4,indxSVSTR)='surface v-mom. stress, scalar, series   '
      vname(5,indxSVSTR)='surface_downward_y_stress               '
      vname(6,indxSVSTR)='lat_v lon_v                             '
      vname(7,indxSVSTR)='                                        '
      return
      end
