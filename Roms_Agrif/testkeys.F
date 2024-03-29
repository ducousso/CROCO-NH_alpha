! $Id: testkeys.F 1458 2014-02-03 15:01:25Z gcambon $
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
#ifdef AGRIF
AGRIFISDEFINED
#endif
#ifdef OPENMP
OPENMPISDEFINED
#endif
#ifdef MPI
MPIISDEFINED
#endif
#if defined OA_COUPLING || defined OW_COUPLING
OACPLISDEFINED
# ifdef OA_MCT
OACPL_MCTISDEFINED
# else
OACPL_OA3ISDEFINED
# endif
#endif
#ifdef XIOS
XIOSISDEFINED
#endif
#ifdef NBQ
NBQISDEFINED
#endif


