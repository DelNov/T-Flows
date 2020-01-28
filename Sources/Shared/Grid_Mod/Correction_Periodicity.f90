!==============================================================================!
  subroutine Grid_Mod_Correction_Periodicity(grid, s, corr_x, corr_y, corr_z)
!------------------------------------------------------------------------------!
!           Find corrections for periodic faces centorid location              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: s
  real            :: corr_x, corr_y, corr_z   ! corrections
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2                           ! counters
!==============================================================================!

  ! Initialize
  corr_x = 0.0; corr_y = 0.0; corr_z = 0.0

  c1 = grid % faces_c(1,s)
  c2 = grid % faces_c(2,s)

  ! In x direction
  if( (abs(grid % per_x) > NANO) .and.                             &
      Math_Mod_Approx_Real(   abs(grid % dx(s))                    &
                            + abs(grid % xc(c1) - grid % xc(c2)),  &
                              grid % per_x) ) then
      corr_x = grid % per_x
  end if

  ! In y direction
  if( (abs(grid % per_y) > NANO) .and.                             &
      Math_Mod_Approx_Real(   abs(grid % dy(s))                    &
                            + abs(grid % yc(c1) - grid % yc(c2)),  &
                              grid % per_y) ) then
      corr_y = grid % per_y
  end if

  ! In z direction
  if( (abs(grid % per_z) > NANO) .and.                             &
      Math_Mod_Approx_Real(   abs(grid % dz(s))                    &
                            + abs(grid % zc(c1) - grid % zc(c2)),  &
                              grid % per_z) ) then
      corr_z = grid % per_z
  end if

  end subroutine
