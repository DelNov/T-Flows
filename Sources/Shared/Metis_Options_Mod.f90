#include "../Shared/Assert.h90"

!==============================================================================!
  module Metis_Options_Mod
!------------------------------------------------------------------------------!
!   Holds parameters for passing options to METIS library.                     !
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
  !==============================================================================!
    subroutine Create_Metis(Metis, face_to_cell, n_parts)
  !------------------------------------------------------------------------------!
    implicit none

    class(Metis_Type)                                :: Metis
    integer, allocatable, intent(in), dimension(:,:) :: face_to_cell
    integer, optional,    intent(in)                 :: n_parts

    integer :: s, i, v, v1, v2

    Assert(n_parts > 0)
  !==============================================================================!

    !------------------------------------------------------------!
    !   Number of vertices and number of edges for first level   !
    !------------------------------------------------------------!
    Metis % n_verts = 0
    Metis % n_edges = 0
    do s = 1, size(face_to_cell, 2)
      if(face_to_cell(2, s) > 0) then
        Metis % n_edges = Metis % n_edges + 1
        Metis % n_verts = max(Metis % n_verts, face_to_cell(2, s))
      end if
    end do
    Assert(Metis % n_verts > 0)
    Assert(Metis % n_edges > 0)

    !---------------------------------------------------------------!
    !   Once n_verts(1) and n_edegs(1) are known, allocate memory   !
    !---------------------------------------------------------------!
    allocate(Metis % edges_v  ( 2,Metis % n_edges));  Metis % edges_v  (:,:) = 0
    allocate(Metis % edges_c  (   Metis % n_edges));  Metis % edges_c  (:)   = 0
    allocate(Metis % star_size(   Metis % n_verts));  Metis % star_size(:)   = 0
    allocate(Metis % star     (24,Metis % n_verts));  Metis % star     (:,:) = 0

    !----------------------------!
    !   Form edge connectivity   !
    !----------------------------!
    i = 0
    do s = 1, size(face_to_cell, 2)
      if(face_to_cell(2, s) > 0) then
        i = i + 1
        Metis % edges_v(1:2, i) = face_to_cell(1:2, s)
        Assert(Metis % edges_v(1,i) > 0)
        Assert(Metis % edges_v(2,i) > 0)
      end if
    end do

    !------------------------------!
    !   Form stars at this level   !
    !------------------------------!
    do s = 1, Metis % n_edges
      v1 = Metis % edges_v(1, s)
      v2 = Metis % edges_v(2, s)

      Metis % star_size(v1) = Metis % star_size(v1) + 1
      Metis % star_size(v2) = Metis % star_size(v2) + 1

      Metis % star(Metis % star_size(v1), v1) = v2
      Metis % star(Metis % star_size(v2), v2) = v1
    end do

    !---------------------------------------------------------------!
    !   Fill-up the structures needed to call METIS (row and col)   !
    !---------------------------------------------------------------!

    ! Fill up the rows
    allocate(Metis % row(Metis % n_verts + 1))
    Metis % row(1) = 0
    do v = 1, Metis % n_verts
      Metis % row(v+1) = Metis % row(v) + Metis % star_size(v)
    end do

    ! Fill up columns
    allocate(Metis % col(Metis % row(Metis % n_verts + 1)))
    do v = 1, Metis % n_verts
      do i = 1, Metis % star_size(v)
        ! In the line below -1 is used because METIS works from 0
        Metis % col(Metis % row(v) + i) = Metis % star(i, v) - 1
      end do
    end do

    !-------------------------------!
    !   Deal with constrains next   !
    !-------------------------------!
    Metis % n_constrains =  1
    allocate(Metis % imbalance(Metis % n_constrains))
    allocate(Metis % vert_weights(Metis % n_verts*Metis % n_constrains))
    allocate(Metis % edge_weights(Metis % row(Metis % n_verts+1)))
    allocate(Metis % vert_data(Metis % n_verts))
    allocate(Metis % part_weight(n_parts*Metis % n_constrains))
    Metis % imbalance(:)    = 1.001
    Metis % vert_weights(:) = 1
    Metis % edge_weights(:) = 1
    Metis % vert_data(:)    = 1
    Metis % part_weight(:)  = 1.0 / real(n_parts)

    end subroutine

  end module
