! $Id: timers_roms.F 1458 2014-02-03 15:01:25Z gcambon $
!
!======================================================================
! ROMS_AGRIF is a branch of ROMS developped at IRD and INRIA, in France
! The two other branches from UCLA (Shchepetkin et al) 
! and Rutgers University (Arango et al) are under MIT/X style license.
! ROMS_AGRIF specific routines (nesting) are under CeCILL-C license.
! 
! ROMS_AGRIF website : http://www.romsagrif.org
!======================================================================
!
#include "cppdefs.h"
      subroutine start_timers()
      implicit none
#if defined SGI || defined SUN
      real*4 etime
      integer  trd, omp_get_thread_num, getpid
#include "param.h"
#include "scalars.h"
!
! Initialize timing for all threads.
!
      trd=omp_get_thread_num()
      if (proc(1,trd).eq.0) then
        proc(1,trd)=1
        proc(0,trd)=getpid()
        timer(0,trd)=etime(timer(1,trd))
C$OMP CRITICAL (start_timers_cr_rgn)
# ifdef MPI
        write(stdout,'(A8,I3,A7,I3,A5,i10,A12)') 'Process ',
     &  mynode,' thread',trd,' (pid=',proc(0,trd),') is active.'
# else
        write(stdout,'(A8,I3,A5,i10,A12)') 'Thread #',trd,
     &                       ' (pid=',proc(0,trd),') is active.'
# endif
C$OMP END CRITICAL (start_timers_cr_rgn)
      endif
#endif
      return
      end



      subroutine stop_timers()
      implicit none
#if defined SGI || defined SUN
      real*4 etime
      integer  trd, omp_get_thread_num
#include "param.h"
#include "scalars.h"
!
! Finalize timings for all threads.
!
      trd=omp_get_thread_num()
      if (proc(1,trd).eq.1) then
        proc(1,trd)=0
        timer(0,trd)=etime(timer(1,trd))-timer(0,trd)
C$OMP CRITICAL (stop_timers_cr_rgn)
# ifdef MPI
         write(stdout,'(A8,I3,A7,I3,A5,F12.3,A5,F10.3,A5,F12.5)')
     &     'Process ',mynode,' thread',trd,' CPU:',timer(1,trd),
     &             ' sys:',timer(2,trd),' net:',timer(0,trd)
# else
         write(stdout,'(A8,I3,A5,F12.3,A5,F10.3,A5,F12.5)')
     &                      'Thread #',trd,' CPU:',timer(1,trd),
     &             ' sys:',timer(2,trd),' net:',timer(0,trd)
# endif
         trd_count=trd_count+1
         timer(0,NPP)=timer(0,NPP)+timer(0,trd)
         timer(1,NPP)=timer(1,NPP)+timer(1,trd)
         timer(2,NPP)=timer(2,NPP)+timer(2,trd)
         if (trd_count.eq.numthreads) then
           trd_count=0
           write(stdout,'(A6,8x,F14.3,3x,F12.3,1x,F14.3)') 'TOTAL:',
     &               timer(1,NPP),timer(2,NPP),timer(0,NPP)
         endif 
C$OMP END CRITICAL (stop_timers_cr_rgn)
      endif
#endif
      return
      end

