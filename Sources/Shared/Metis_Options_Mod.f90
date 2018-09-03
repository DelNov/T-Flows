!==============================================================================!
  module Metis_Options_Mod
!------------------------------------------------------------------------------!
!   Holds parameters for passing options to METIS library.                     !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  integer, parameter :: METIS_OPTION_PTYPE     =  1
  integer, parameter :: METIS_OPTION_OBJTYPE   =  2
  integer, parameter :: METIS_OPTION_CTYPE     =  3
  integer, parameter :: METIS_OPTION_IPTYPE    =  4
  integer, parameter :: METIS_OPTION_RTYPE     =  5
  integer, parameter :: METIS_OPTION_DBGLVL    =  6
  integer, parameter :: METIS_OPTION_NITER     =  7
  integer, parameter :: METIS_OPTION_NCUTS     =  8
  integer, parameter :: METIS_OPTION_SEED      =  9
  integer, parameter :: METIS_OPTION_NO2HOP    = 10
  integer, parameter :: METIS_OPTION_MINCONN   = 11
  integer, parameter :: METIS_OPTION_CONTIG    = 12
  integer, parameter :: METIS_OPTION_COMPRESS  = 13
  integer, parameter :: METIS_OPTION_CCORDER   = 14
  integer, parameter :: METIS_OPTION_PFACTOR   = 15
  integer, parameter :: METIS_OPTION_NSEPS     = 16
  integer, parameter :: METIS_OPTION_UFACTOR   = 17
  integer, parameter :: METIS_OPTION_NUMBERING = 18
  integer, parameter :: METIS_OPTION_HELP      = 19
  integer, parameter :: METIS_OPTION_TPWGTS    = 20
  integer, parameter :: METIS_OPTION_NCOMMON   = 21
  integer, parameter :: METIS_OPTION_NOOUTPUT  = 22
  integer, parameter :: METIS_OPTION_BALANCE   = 23
  integer, parameter :: METIS_OPTION_GTYPE     = 24
  integer, parameter :: METIS_OPTION_UBVEC     = 25

  end module
