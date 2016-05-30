      subroutine MessPass2D_tile (Istr,Iend,Jstr,Jend, A)
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      include 'mpif.h'
      integer*4 Npts,ipts,jpts
      parameter (Npts=2)
      real A(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      integer*4 Istr,Iend,Jstr,Jend, i,j, isize,jsize,ksize,  iter,
     &        req(8), status(MPI_STATUS_SIZE,8), ierr, mdii,mdjj
      real buf_snd4(Npts*Npts),     buf_snd2(Npts*Npts),
     &     buf_rev4(Npts*Npts),     buf_rev2(Npts*Npts),
     &     buf_snd1(Npts*Npts),     buf_snd3(Npts*Npts),
     &     buf_rev1(Npts*Npts),     buf_rev3(Npts*Npts)
      integer*4 sub_X,size_X, sub_E,size_E
      parameter (sub_X=Lm,  size_X=7+Npts*sub_X,
     &           sub_E=Mm,  size_E=7+Npts*sub_E)
      real ibuf_sndN(0:size_X), ibuf_revN(0:size_X),
     &     ibuf_sndS(0:size_X), ibuf_revS(0:size_X),
     &     jbuf_sndW(0:size_E), jbuf_sndE(0:size_E),
     &     jbuf_revW(0:size_E), jbuf_revE(0:size_E)
      integer*4 imin,imax,ishft, jmin,jmax,jshft
      if (ii.eq.0 .and. Istr.eq.1) then
        imin=Istr-1
      else
        imin=Istr
      endif
      if (ii.eq.NP_XI-1 .and. Iend.eq.Lmmpi) then
        imax=Iend+1
      else
        imax=Iend
      endif
      ishft=imax-imin+1
      if (jj.eq.0 .and. Jstr.eq.1) then
        jmin=Jstr-1
      else
        jmin=Jstr
      endif
      if (jj.eq.NP_ETA-1 .and. Jend.eq.Mmmpi) then
        jmax=Jend+1
      else
        jmax=Jend
      endif
      jshft=jmax-jmin+1
      ksize=Npts*Npts
      isize=Npts*ishft
      jsize=Npts*jshft
      do iter=0,1
        mdii=mod(ii+iter,2)
        mdjj=mod(jj+iter,2)
        if (mdii.eq.0) then
          if (WEST_INTER) then
            do j=jmin,jmax
              do ipts=1,Npts
                jbuf_sndW(j-jmin+(ipts-1)*jshft)=A(ipts,j)
              enddo
            enddo
            call MPI_Irecv (jbuf_revW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 2, MPI_COMM_WORLD, req(1), ierr)
            call MPI_Send  (jbuf_sndW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 1, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (EAST_INTER) then
            do j=jmin,jmax
              do ipts=1,Npts
                jbuf_sndE(j-jmin+(ipts-1)*jshft)=A(Lmmpi-Npts+ipts,j)
              enddo
            enddo
            call MPI_Irecv (jbuf_revE, jsize, MPI_DOUBLE_PRECISION,
     &                         p_E, 1, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Send  (jbuf_sndE, jsize, MPI_DOUBLE_PRECISION,
     &                         p_E, 2, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdjj.eq.0) then
          if (SOUTH_INTER) then
            do i=imin,imax
              do jpts=1,Npts
                ibuf_sndS(i-imin+(jpts-1)*ishft)=A(i,jpts)
              enddo
            enddo
            call MPI_Irecv (ibuf_revS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 4, MPI_COMM_WORLD, req(3), ierr)
            call MPI_Send  (ibuf_sndS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 3, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER) then
            do i=imin,imax
              do jpts=1,Npts
                ibuf_sndN(i-imin+(jpts-1)*ishft)=A(i,Mmmpi-Npts+jpts)
              enddo
            enddo
            call MPI_Irecv (ibuf_revN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 3, MPI_COMM_WORLD, req(4), ierr)
            call MPI_Send  (ibuf_sndN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 4, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdii.eq.0) then
          if (SOUTH_INTER .and. WEST_INTER) then
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd1(ipts+Npts*(jpts-1))=A(ipts,jpts)
              enddo
            enddo
            call MPI_Irecv (buf_rev1,ksize, MPI_DOUBLE_PRECISION,  p_SW,
     &                               6, MPI_COMM_WORLD, req(5),ierr)
            call MPI_Send  (buf_snd1,ksize, MPI_DOUBLE_PRECISION,  p_SW,
     &                               5, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (NORTH_INTER .and. EAST_INTER) then
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd2(ipts+Npts*(jpts-1))=
     &                                A(Lmmpi-Npts+ipts,Mmmpi-Npts+jpts)
              enddo
            enddo
            call MPI_Irecv (buf_rev2,ksize, MPI_DOUBLE_PRECISION,  p_NE,
     &                               5, MPI_COMM_WORLD, req(6),ierr)
            call MPI_Send  (buf_snd2,ksize, MPI_DOUBLE_PRECISION,  p_NE,
     &                               6, MPI_COMM_WORLD,        ierr)
          endif
        endif
        if (mdii.eq.1) then
          if (SOUTH_INTER .and. EAST_INTER) then
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd3(ipts+Npts*(jpts-1))=
     &                                A(Lmmpi-Npts+ipts,jpts)
              enddo
            enddo
            call MPI_Irecv (buf_rev3,ksize, MPI_DOUBLE_PRECISION,  p_SE,
     &                               8, MPI_COMM_WORLD, req(7), ierr)
            call MPI_Send  (buf_snd3,ksize, MPI_DOUBLE_PRECISION,  p_SE,
     &                               7, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER .and. WEST_INTER) then
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd4(ipts+Npts*(jpts-1))=
     &                                A(ipts,Mmmpi-Npts+jpts)
              enddo
            enddo
            call MPI_Irecv (buf_rev4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                               7, MPI_COMM_WORLD, req(8), ierr)
            call MPI_Send  (buf_snd4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                               8, MPI_COMM_WORLD,         ierr)
          endif
        endif
      enddo
      if (WEST_INTER) then
        call MPI_Wait (req(1),status(1,1),ierr)
        do j=jmin,jmax
          do ipts=1,Npts
           A(ipts-Npts,j)=jbuf_revW(j-jmin+(ipts-1)*jshft)
          enddo
        enddo
      endif
      if (EAST_INTER) then
        call MPI_Wait (req(2),status(1,2),ierr)
        do j=jmin,jmax
          do ipts=1,Npts
           A(Lmmpi+ipts,j)=jbuf_revE(j-jmin+(ipts-1)*jshft)
          enddo
        enddo
      endif
      if (SOUTH_INTER) then
        call MPI_Wait (req(3),status(1,3),ierr)
        do i=imin,imax
          do jpts=1,Npts
           A(i,jpts-Npts)=ibuf_revS(i-imin+(jpts-1)*ishft)
          enddo
        enddo
      endif
      if (NORTH_INTER) then
        call MPI_Wait (req(4),status(1,4),ierr)
        do i=imin,imax
          do jpts=1,Npts
            A(i,Mmmpi+jpts)=ibuf_revN(i-imin+(jpts-1)*ishft)
          enddo
        enddo
      endif
      if (SOUTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(5),status(1,5),ierr)
        do jpts=1,Npts
          do ipts=1,Npts
            A(ipts-Npts,jpts-Npts)=buf_rev1(ipts+Npts*(jpts-1))
          enddo
        enddo
      endif
      if (NORTH_INTER .and. EAST_INTER) then
        call MPI_Wait (req(6),status(1,6),ierr)
        do jpts=1,Npts
          do ipts=1,Npts
            A(Lmmpi+ipts,Mmmpi+jpts)=buf_rev2(ipts+Npts*(jpts-1))
          enddo
        enddo
      endif
      if (SOUTH_INTER .and. EAST_INTER) then
        call MPI_Wait (req(7),status(1,7),ierr)
        do jpts=1,Npts
          do ipts=1,Npts
            A(Lmmpi+ipts,jpts-Npts)=buf_rev3(ipts+Npts*(jpts-1))
          enddo
        enddo
      endif
      if (NORTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(8),status(1,8),ierr)
        do jpts=1,Npts
          do ipts=1,Npts
            A(ipts-Npts,Mmmpi+jpts)=buf_rev4(ipts+Npts*(jpts-1))
          enddo
        enddo
      endif
      return
      end
      subroutine MessPass2D_3pts_tile (Istr,Iend,Jstr,Jend, A)
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      include 'mpif.h'
      integer*4 Npts,ipts,jpts
      parameter (Npts=3)
      real A(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      integer*4 Istr,Iend,Jstr,Jend, i,j, isize,jsize,ksize,  iter,
     &        req(8), status(MPI_STATUS_SIZE,8), ierr, mdii,mdjj
      real buf_snd4(Npts*Npts),     buf_snd2(Npts*Npts),
     &     buf_rev4(Npts*Npts),     buf_rev2(Npts*Npts),
     &     buf_snd1(Npts*Npts),     buf_snd3(Npts*Npts),
     &     buf_rev1(Npts*Npts),     buf_rev3(Npts*Npts)
      integer*4 sub_X,size_X, sub_E,size_E
      parameter (sub_X=Lm,  size_X=7+Npts*sub_X,
     &           sub_E=Mm,  size_E=7+Npts*sub_E)
      real ibuf_sndN(0:size_X), ibuf_revN(0:size_X),
     &     ibuf_sndS(0:size_X), ibuf_revS(0:size_X),
     &     jbuf_sndW(0:size_E), jbuf_sndE(0:size_E),
     &     jbuf_revW(0:size_E), jbuf_revE(0:size_E)
      integer*4 imin,imax,ishft, jmin,jmax,jshft
      if (ii.eq.0 .and. Istr.eq.1) then
        imin=Istr-1
      else
        imin=Istr
      endif
      if (ii.eq.NP_XI-1 .and. Iend.eq.Lmmpi) then
        imax=Iend+1
      else
        imax=Iend
      endif
      ishft=imax-imin+1
      if (jj.eq.0 .and. Jstr.eq.1) then
        jmin=Jstr-1
      else
        jmin=Jstr
      endif
      if (jj.eq.NP_ETA-1 .and. Jend.eq.Mmmpi) then
        jmax=Jend+1
      else
        jmax=Jend
      endif
      jshft=jmax-jmin+1
      ksize=Npts*Npts
      isize=Npts*ishft
      jsize=Npts*jshft
      do iter=0,1
        mdii=mod(ii+iter,2)
        mdjj=mod(jj+iter,2)
        if (mdii.eq.0) then
          if (WEST_INTER) then
            do j=jmin,jmax
              do ipts=1,Npts
                jbuf_sndW(j-jmin+(ipts-1)*jshft)=A(ipts,j)
              enddo
            enddo
            call MPI_Irecv (jbuf_revW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 2, MPI_COMM_WORLD, req(1), ierr)
            call MPI_Send  (jbuf_sndW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, 1, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (EAST_INTER) then
            do j=jmin,jmax
              do ipts=1,Npts
                jbuf_sndE(j-jmin+(ipts-1)*jshft)=A(Lmmpi-Npts+ipts,j)
              enddo
            enddo
            call MPI_Irecv (jbuf_revE, jsize, MPI_DOUBLE_PRECISION,
     &                         p_E, 1, MPI_COMM_WORLD, req(2), ierr)
            call MPI_Send  (jbuf_sndE, jsize, MPI_DOUBLE_PRECISION,
     &                         p_E, 2, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdjj.eq.0) then
          if (SOUTH_INTER) then
            do i=imin,imax
              do jpts=1,Npts
                ibuf_sndS(i-imin+(jpts-1)*ishft)=A(i,jpts)
              enddo
            enddo
            call MPI_Irecv (ibuf_revS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 4, MPI_COMM_WORLD, req(3), ierr)
            call MPI_Send  (ibuf_sndS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, 3, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER) then
            do i=imin,imax
              do jpts=1,Npts
                ibuf_sndN(i-imin+(jpts-1)*ishft)=A(i,Mmmpi-Npts+jpts)
              enddo
            enddo
            call MPI_Irecv (ibuf_revN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 3, MPI_COMM_WORLD, req(4), ierr)
            call MPI_Send  (ibuf_sndN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, 4, MPI_COMM_WORLD,         ierr)
          endif
        endif
        if (mdii.eq.0) then
          if (SOUTH_INTER .and. WEST_INTER) then
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd1(ipts+Npts*(jpts-1))=A(ipts,jpts)
              enddo
            enddo
            call MPI_Irecv (buf_rev1,ksize, MPI_DOUBLE_PRECISION,  p_SW,
     &                               6, MPI_COMM_WORLD, req(5),ierr)
            call MPI_Send  (buf_snd1,ksize, MPI_DOUBLE_PRECISION,  p_SW,
     &                               5, MPI_COMM_WORLD,        ierr)
          endif
        else
          if (NORTH_INTER .and. EAST_INTER) then
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd2(ipts+Npts*(jpts-1))=
     &                                A(Lmmpi-Npts+ipts,Mmmpi-Npts+jpts)
              enddo
            enddo
            call MPI_Irecv (buf_rev2,ksize, MPI_DOUBLE_PRECISION,  p_NE,
     &                               5, MPI_COMM_WORLD, req(6),ierr)
            call MPI_Send  (buf_snd2,ksize, MPI_DOUBLE_PRECISION,  p_NE,
     &                               6, MPI_COMM_WORLD,        ierr)
          endif
        endif
        if (mdii.eq.1) then
          if (SOUTH_INTER .and. EAST_INTER) then
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd3(ipts+Npts*(jpts-1))=
     &                                A(Lmmpi-Npts+ipts,jpts)
              enddo
            enddo
            call MPI_Irecv (buf_rev3,ksize, MPI_DOUBLE_PRECISION,  p_SE,
     &                               8, MPI_COMM_WORLD, req(7), ierr)
            call MPI_Send  (buf_snd3,ksize, MPI_DOUBLE_PRECISION,  p_SE,
     &                               7, MPI_COMM_WORLD,         ierr)
          endif
        else
          if (NORTH_INTER .and. WEST_INTER) then
            do jpts=1,Npts
              do ipts=1,Npts
                buf_snd4(ipts+Npts*(jpts-1))=
     &                                A(ipts,Mmmpi-Npts+jpts)
              enddo
            enddo
            call MPI_Irecv (buf_rev4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                               7, MPI_COMM_WORLD, req(8), ierr)
            call MPI_Send  (buf_snd4, ksize, MPI_DOUBLE_PRECISION, p_NW,
     &                               8, MPI_COMM_WORLD,         ierr)
          endif
        endif
      enddo
      if (WEST_INTER) then
        call MPI_Wait (req(1),status(1,1),ierr)
        do j=jmin,jmax
          do ipts=1,Npts
           A(ipts-Npts,j)=jbuf_revW(j-jmin+(ipts-1)*jshft)
          enddo
        enddo
      endif
      if (EAST_INTER) then
        call MPI_Wait (req(2),status(1,2),ierr)
        do j=jmin,jmax
          do ipts=1,Npts
           A(Lmmpi+ipts,j)=jbuf_revE(j-jmin+(ipts-1)*jshft)
          enddo
        enddo
      endif
      if (SOUTH_INTER) then
        call MPI_Wait (req(3),status(1,3),ierr)
        do i=imin,imax
          do jpts=1,Npts
           A(i,jpts-Npts)=ibuf_revS(i-imin+(jpts-1)*ishft)
          enddo
        enddo
      endif
      if (NORTH_INTER) then
        call MPI_Wait (req(4),status(1,4),ierr)
        do i=imin,imax
          do jpts=1,Npts
            A(i,Mmmpi+jpts)=ibuf_revN(i-imin+(jpts-1)*ishft)
          enddo
        enddo
      endif
      if (SOUTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(5),status(1,5),ierr)
        do jpts=1,Npts
          do ipts=1,Npts
            A(ipts-Npts,jpts-Npts)=buf_rev1(ipts+Npts*(jpts-1))
          enddo
        enddo
      endif
      if (NORTH_INTER .and. EAST_INTER) then
        call MPI_Wait (req(6),status(1,6),ierr)
        do jpts=1,Npts
          do ipts=1,Npts
            A(Lmmpi+ipts,Mmmpi+jpts)=buf_rev2(ipts+Npts*(jpts-1))
          enddo
        enddo
      endif
      if (SOUTH_INTER .and. EAST_INTER) then
        call MPI_Wait (req(7),status(1,7),ierr)
        do jpts=1,Npts
          do ipts=1,Npts
            A(Lmmpi+ipts,jpts-Npts)=buf_rev3(ipts+Npts*(jpts-1))
          enddo
        enddo
      endif
      if (NORTH_INTER .and. WEST_INTER) then
        call MPI_Wait (req(8),status(1,8),ierr)
        do jpts=1,Npts
          do ipts=1,Npts
            A(ipts-Npts,Mmmpi+jpts)=buf_rev4(ipts+Npts*(jpts-1))
          enddo
        enddo
      endif
      return
      end
