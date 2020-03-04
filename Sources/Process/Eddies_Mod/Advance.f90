!==============================================================================!
  subroutine Eddies_Mod_Advance(eddies)
!------------------------------------------------------------------------------!
!   Advance all eddies through the inlet plane                                 !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Eddies_Type), target :: eddies
  real                      :: dt      ! global time step
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  integer                   :: e
!==============================================================================!

  ! Create an alias
  grid => eddies % pnt_grid
  flow => eddies % pnt_flow

  ! Advance all eddies
  do e = 1, eddies % n_eddies
    eddies % eddy(e) % x = eddies % eddy(e) % x  &
                         + eddies % eddy(e) % u * flow % dt
    eddies % eddy(e) % y = eddies % eddy(e) % y  &
                         + eddies % eddy(e) % v * flow % dt
    eddies % eddy(e) % z = eddies % eddy(e) % z  &
                         + eddies % eddy(e) % w * flow % dt
  end do

  ! Check which eddies left
  if(eddies % x_plane < HUGE) then
    do e = 1, eddies % n_eddies
      if( abs(eddies % eddy(e) % x - eddies % x_plane)  &
            > eddies % eddy(e) % length) then
        call Eddies_Mod_Place_Eddy(eddies, e)
      end if
    end do
  else if(eddies % y_plane < HUGE) then
    do e = 1, eddies % n_eddies
      if( abs(eddies % eddy(e) % y - eddies % y_plane)  &
            > eddies % eddy(e) % length) then
        call Eddies_Mod_Place_Eddy(eddies, e)
      end if
    end do
  else if(eddies % z_plane < HUGE) then
    do e = 1, eddies % n_eddies
      if( abs(eddies % eddy(e) % z - eddies % z_plane)  &
            > eddies % eddy(e) % length) then
        call Eddies_Mod_Place_Eddy(eddies, e)
      end if
    end do
  end if

  end subroutine
