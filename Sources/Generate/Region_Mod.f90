!==============================================================================!
  module Region_Mod
!------------------------------------------------------------------------------!
!   This is used to read boundary conditions or materials in "Generator"       !
!------------------------------------------------------------------------------!
  implicit none
!==============================================================================!

  !-----------------!
  !   Region type   !
  !-----------------!
  type Region_Type

    integer :: block   ! to which block it belons

    ! Enclosing logical cooridnates
    integer :: is, ie  ! start and end in "i" direction
    integer :: js, je  ! start and end in "i" direction
    integer :: ks, ke  ! start and end in "i" direction

    ! On which face is it ("IMIN", "IMAC", "JMIN", ...
    character(len= 4) :: face

    ! Stores the name of the boundary condition (or material)
    character(len=80) :: name

  end type

  end module
