!==============================================================================!
  module Save_Grid_Mod
!------------------------------------------------------------------------------!
!   Module for saving results.                                                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Grid_Mod
  use Cgns_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Constants shared among most of the procedures in this module
  integer,           parameter :: VTK_LINE       =  3
  integer,           parameter :: VTK_TRIANGLE   =  5
  integer,           parameter :: VTK_QUAD       =  9
  integer,           parameter :: VTK_TETRA      = 10
  integer,           parameter :: VTK_HEXAHEDRON = 12
  integer,           parameter :: VTK_WEDGE      = 13
  integer,           parameter :: VTK_PYRAMID    = 14
  character(len= 1), parameter :: LF   = char(10)      ! line feed
  character(len= 0), parameter :: IN_0 = ''            ! indentation levels
  character(len= 2), parameter :: IN_1 = '  '
  character(len= 4), parameter :: IN_2 = '    '
  character(len= 6), parameter :: IN_3 = '      '
  character(len= 8), parameter :: IN_4 = '        '
  character(len=10), parameter :: IN_5 = '          '

  contains

  include 'Save_Grid_Mod/Vtu/Save_Vtu_Cells.f90'
  include 'Save_Grid_Mod/Vtu/Save_Vtu_Faces.f90'
  include 'Save_Grid_Mod/Cgns/Save_Cgns_Cells.f90'

  end module 
