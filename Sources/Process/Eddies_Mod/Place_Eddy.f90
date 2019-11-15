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
! if(min_dist < ONE_THIRD * eddies % max_radius) then
!   goto 1
! end if

  ! Assume eddies move in x direction
  call random_number(tmp)
  eddies % eddy(e) % x = eddies % eddy(e) % x             &
                       - eddies % eddy(e) % length * 0.5  &
                       + eddies % eddy(e) % length * tmp

!DEBUG EDDIES % EDDY(E) % Y = 0.5
!DEBUG EDDIES % EDDY(E) % Z = 0.5
!DEBUG EDDIES % EDDY(E) % RADIUS = 0.45
!DEBUG EDDIES % EDDY(E) % LENGTH = 3.0

  end subroutine
