!==============================================================================!
  module Range_Mod
!------------------------------------------------------------------------------!
!>  This module is used for defining and handling boundary conditions or
!>  materials for the computational grids created with Generate.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Range type   !
  !----------------!
  !> Central to the Range_Mod. It encapsulates the necessary
  !> information to define a range within a computational grid.
  type Range_Type

    integer :: block   !! for which block is the range defined

    ! Enclosing logical cooridnates
    integer :: is, ie  !! start and end in "i" direction
    integer :: js, je  !! start and end in "j" direction
    integer :: ks, ke  !! start and end in "k" direction

    ! On which face is it ("IMIN", "IMAC", "JMIN", ...
    character(SL) :: face  !! string indicating on which face of the block the
                           !! range is located (IMIN, IMAX, JMIN, ... ZMAX)

    ! Stores the name of the boundary condition (or material)
    character(SL) :: name  !! string to store the name of the region
                           !! (boundary condition or material)
  end type

  end module
