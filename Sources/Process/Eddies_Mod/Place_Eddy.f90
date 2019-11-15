!==============================================================================!
  subroutine Eddies_Mod_Place_Eddy(eddies, e)
!------------------------------------------------------------------------------!
!   Place one eddy randomly                                                    !
!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Eddies_Type), target :: eddies
  integer                   :: e
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: ri, rc, oe
  real                     :: tmp, min_dist
!==============================================================================!

  ! Create an alias
  grid => eddies % pnt_grid

  !-----------------------!
  !   Set eddy's radius   !
  !-----------------------!
  call random_number(tmp)
  eddies % eddy(e) % radius = 0.25 * eddies % max_radius  &
                            + 0.75 * eddies % max_radius * tmp
  eddies % eddy(e) % length = eddies % eddy(e) % radius * 6.0

  !-------------------------!
  !   Set eddy's rotation   !
  !-------------------------!
  call random_number(tmp)
  eddies % eddy(e) % sense = 1.0
  if(tmp < 0.5) eddies % eddy(e) % sense = -1.0

  !-----------------------------------------------!
  !   Position eddy in the random boundary cell   !
  !-----------------------------------------------!
1 continue
  call random_number(tmp)
  rc = int(tmp * real(eddies % n_bnd_cells_glo))  ! random index
  eddies % eddy(e) % x = eddies % bnd_xc(rc)
  eddies % eddy(e) % y = eddies % bnd_yc(rc)
  eddies % eddy(e) % z = eddies % bnd_zc(rc)
  eddies % eddy(e) % u = eddies % bnd_u (rc)
  eddies % eddy(e) % v = eddies % bnd_v (rc)
  eddies % eddy(e) % w = eddies % bnd_w (rc)

  ! Don't allow eddies bigger than wall distance
  if(eddies % eddy(e) % radius > eddies % bnd_wd(rc)) goto 1

  ! Check that it is not too close to other eddies
  min_dist = HUGE
  do oe = 1, eddies % n_eddies
    if(oe .ne. e) then
      min_dist = min(min_dist, Math_Mod_Distance(eddies % eddy(e)  % x,  &
                                                 eddies % eddy(e)  % y,  &
                                                 eddies % eddy(e)  % z,  &
                                                 eddies % eddy(oe) % x,  &
                                                 eddies % eddy(oe) % y,  &
                                                 eddies % eddy(oe) % z))
    end if
  end do
  if(min_dist < TWO_THIRDS * eddies % max_radius) then
    goto 1
  end if

  ! Assume eddies move in x direction
  call random_number(tmp)
  if(eddies % x_plane < HUGE) then
    eddies % eddy(e) % x = eddies % eddy(e) % x             &
                         - eddies % eddy(e) % length * 0.5  &
                         + eddies % eddy(e) % length * tmp
  else if(eddies % y_plane < HUGE) then
    eddies % eddy(e) % y = eddies % eddy(e) % y             &
                         - eddies % eddy(e) % length * 0.5  &
                         + eddies % eddy(e) % length * tmp
  else if(eddies % z_plane < HUGE) then
    eddies % eddy(e) % z = eddies % eddy(e) % z             &
                         - eddies % eddy(e) % length * 0.5  &
                         + eddies % eddy(e) % length * tmp
  end if

  end subroutine
