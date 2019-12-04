!==============================================================================!
  module Bnd_Cond_Mod
!------------------------------------------------------------------------------!
!   This is used to store boundary conditions within a Grid_Type               !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !---------------------------------------------------------!
  !   Constants for identification of boundary conditions   !
  !---------------------------------------------------------!
  integer, parameter :: INFLOW   = 10007
  integer, parameter :: WALL     = 10009
  integer, parameter :: OUTFLOW  = 10037
  integer, parameter :: SYMMETRY = 10039
  integer, parameter :: CONVECT  = 10061
  integer, parameter :: WALLFL   = 10067
  integer, parameter :: PRESSURE = 10069
  integer, parameter :: PERIODIC = 10079

  !-------------------!
  !   Bnd_Cond type   !
  !-------------------!
  type Bnd_Cond_Type

    ! Name of the boundary conditions specified in grid generation
    ! It ranges through number of boundary conditions (aka colors)
    character(len=80), allocatable :: name(:)

    ! Boundary types, ranging through all colors
    integer, allocatable :: type(:)

    ! Boundary condition color ranging through boundary cells.
    ! Values start from one, zero is internal cell
    ! (Follows nomenclature from "../Shared/Comm_Mod_Par.f90")
    integer, allocatable :: color(:)
    integer, allocatable :: color_s_cell(:)  ! start bnd cell for color
    integer, allocatable :: color_e_cell(:)  ! end bnd cell for color

  end type

  end module
