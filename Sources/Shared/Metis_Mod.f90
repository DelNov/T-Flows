#include "../Shared/Assert.h90"

!==============================================================================!
  module Metis_Mod
!------------------------------------------------------------------------------!
!   Holds parameters and a procedure for interaction with METIS library.       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Assert_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Metis type   !
  !----------------!
  type Metis_Type

    ! Data members to prepare a call to METIS
    integer, allocatable :: star     (:,:),  &
                            star_size  (:),  &
                            edges_v  (:,:),  &
                            edges_c    (:)

    ! Data members to make an actual call to METIS
    integer                :: n_verts, n_edges
    integer, dimension(41) :: options              ! options passed to METIS
    integer                :: n_constrains,     &  ! number of constraints
                            return_val             ! return value from METIS
    integer, allocatable :: row(:),             &  ! incidency matrix in ...
                            col(:),             &  ! compresser row format
                            vert_weights(:),    &  ! weights of vertices
                            edge_weights(:),    &  ! weights of edges
                            vert_data(:)           ! amount of data for vertices
    real, allocatable    :: part_weight(:),     &
                            imbalance(:)           ! allowed imbalance

    contains
      procedure :: Create_Metis

  end type

  type(Metis_Type) :: Metis

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

  contains
#   include "Metis_Mod/Create_Metis.f90"

  end module
