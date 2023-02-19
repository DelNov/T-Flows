!==============================================================================!
  module Range_Mod
!------------------------------------------------------------------------------!
!   This is used to read boundary conditions or materials in "Generator"       !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !----------------!
  !   Range type   !
  !----------------!
  type Range_Type

    integer :: block   ! to which block it belons

    ! Enclosing logical cooridnates
    integer :: is, ie  ! start and end in "i" direction
    integer :: js, je  ! start and end in "i" direction
    integer :: ks, ke  ! start and end in "i" direction

    ! On which face is it ("IMIN", "IMAC", "JMIN", ...
    character(SL) :: face

    ! Stores the name of the boundary condition (or material)
    character(SL) :: name

  end type

  end module
