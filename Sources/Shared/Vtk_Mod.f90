!==============================================================================!
  module Vtk_Mod
!------------------------------------------------------------------------------!
!   This is used to store parameters associated with VTK file format           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !------------------------------------!
  !   Constants defining cell shapes   !
  !------------------------------------!
  integer, parameter :: VTK_LINE       =  3
  integer, parameter :: VTK_TRIANGLE   =  5
  integer, parameter :: VTK_POLYGON    =  7
  integer, parameter :: VTK_QUAD       =  9
  integer, parameter :: VTK_TETRA      = 10
  integer, parameter :: VTK_HEXAHEDRON = 12
  integer, parameter :: VTK_WEDGE      = 13
  integer, parameter :: VTK_PYRAMID    = 14
  integer, parameter :: VTK_POLYHEDRON = 42

  !-------------------------------------------------------!
  !   Constants for text formatting; indentation levels   !
  !-------------------------------------------------------!
  character(len= 1), parameter :: LF   = char(10)      ! line feed
  character(len= 0), parameter :: IN_0 = ''            ! indentation levels
  character(len= 2), parameter :: IN_1 = '  '
  character(len= 4), parameter :: IN_2 = '    '
  character(len= 6), parameter :: IN_3 = '      '
  character(len= 8), parameter :: IN_4 = '        '
  character(len=10), parameter :: IN_5 = '          '

  character(len=7) :: intp   = '"IntXX"'
  character(len=9) :: floatp = '"FloatXX"'

  contains

    include 'Vtk_Mod/Set_Precision.f90'

  end module
