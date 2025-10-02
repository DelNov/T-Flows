!==============================================================================!
!   Introduce new types to be used with User_Mod                               !
!==============================================================================!

  ! Indices to store difference species
  integer, parameter :: PB_G   = 1  ! Pb(g)
  integer, parameter :: PBI_G  = 2  ! PbI(g)
  integer, parameter :: I_G    = 3  ! I(g)
  integer, parameter :: PB_S   = 4  ! Pb(s)
  integer, parameter :: PBI_S  = 5  ! PbI(s)
  integer, parameter :: I_S    = 6  ! I(s)

  ! Reaction rates
  real, allocatable :: rate(:,:)

