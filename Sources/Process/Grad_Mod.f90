!==============================================================================!
  module Grad_Mod
!------------------------------------------------------------------------------!
!   Module used for calculation of gradients.                                  !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  ! Gradient matrices for all cells
  real, allocatable :: g(:,:)

  ! Cells which are bad for calculation of gradients
  logical, allocatable :: bad_cells(:)

  contains

  include 'Grad_Mod/Allocate.f90'
  include 'Grad_Mod/Array.f90'
  include 'Grad_Mod/Calculate_Matrix.f90'
  include 'Grad_Mod/Correct_Bad_Cells.f90'
  include 'Grad_Mod/Find_Bad_Cells.f90'
  include 'Grad_Mod/Pressure.f90'
  include 'Grad_Mod/Component.f90'
  include 'Grad_Mod/Variable.f90'

  end module 
