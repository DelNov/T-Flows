!==============================================================================!
  subroutine Place_Eddy(Eddies, e)
!------------------------------------------------------------------------------!
!>  This subroutine is designed to randomly position a single eddy within the
!>  specified boundary in a computational domain used for scale-resolving
!>  simulations.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Eddy characteristics:                                                    !
!     - The subroutine randomly determines the radius and length of the eddy.  !
!       The radius is a fraction of the maximum radius specified in the eddies !
!       structure, and the length is six times the radius.                     !
!     - The sense of rotation (clockwise or counterclockwise) is also randomly !
!       assigned.                                                              !
!   * Random positioning of eddy:                                              !
!     - The subroutine selects a random boundary cell and sets the eddy's      !
!       position (x, y, z) and velocity components (u, v, w) based on the      !
!       boundary cell's properties.                                            !
!     - It ensures that the eddy's radius does not exceed the wall distance    !
!       of the selected boundary cell. If it does, another cell is chosen.     !
!   * Distance check:                                                          !
!     - The subroutine checks that the newly positioned eddy is not too        !
!       close to other existing eddies. If it's within a certain threshold     !
!       distance, a new position is selected.                                  !
!   * Placement relative to inlet plane:                                       !
!     - The eddy is placed with a certain shift relative to the inlet plane    !
!       (x-plane, y-plane, or z-plane). This is done to ensure the eddy enters !
!       the computational domain smoothly.                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Eddies_Type), target :: Eddies  !! parent class; collection of eddies
  integer,        intent(in) :: e       !! index of an eddy
!-----------------------------------[Locals]-----------------------------------!
  integer :: rc, oe
  real    :: tmp, min_dist
!==============================================================================!

  !-----------------------!
  !   Set eddy's radius   !
  !-----------------------!
  call random_number(tmp)
  Eddies % eddy(e) % radius = 0.25 * Eddies % max_radius  &
                            + 0.75 * Eddies % max_radius * tmp
  Eddies % eddy(e) % length = Eddies % eddy(e) % radius * 6.0

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
  rc = 1 + floor(tmp * real(Eddies % n_bnd_cells_glo))  ! random index
  Eddies % eddy(e) % x = Eddies % bnd_xc(rc)
  Eddies % eddy(e) % y = Eddies % bnd_yc(rc)
  Eddies % eddy(e) % z = Eddies % bnd_zc(rc)
  Eddies % eddy(e) % u = Eddies % bnd_u (rc)
  Eddies % eddy(e) % v = Eddies % bnd_v (rc)
  Eddies % eddy(e) % w = Eddies % bnd_w (rc)

  ! Don't allow eddies bigger than wall distance
  if(Eddies % eddy(e) % radius > Eddies % bnd_wd(rc)) goto 1

  ! Check that it is not too close to other eddies
  min_dist = HUGE
  do oe = 1, Eddies % n_eddies
    if(oe .ne. e) then
      min_dist = min(min_dist, Math % Distance(Eddies % eddy(e)  % x,  &
                                               Eddies % eddy(e)  % y,  &
                                               Eddies % eddy(e)  % z,  &
                                               Eddies % eddy(oe) % x,  &
                                               Eddies % eddy(oe) % y,  &
                                               Eddies % eddy(oe) % z))
    end if
  end do
  if(min_dist < TWO_THIRDS * Eddies % max_radius) then
    goto 1
  end if

  ! Place eddies with a certain shift to the inlet plane
  call random_number(tmp)
  if(Eddies % x_plane < HUGE) then
    Eddies % eddy(e) % x = Eddies % eddy(e) % x                 &
                         + sign(1.0, Eddies % eddy(e) % u)      &
                         * (- Eddies % eddy(e) % length * 2.0   &
                            + Eddies % eddy(e) % length * tmp)
  else if(Eddies % y_plane < HUGE) then
    Eddies % eddy(e) % y = Eddies % eddy(e) % y                 &
                         + sign(1.0, Eddies % eddy(e) % v)      &
                         * (- Eddies % eddy(e) % length * 2.0   &
                            + Eddies % eddy(e) % length * tmp)
  else if(Eddies % z_plane < HUGE) then
    Eddies % eddy(e) % z = Eddies % eddy(e) % z                 &
                         + sign(1.0, Eddies % eddy(e) % w)      &
                         * (- Eddies % eddy(e) % length * 2.0   &
                            + Eddies % eddy(e) % length * tmp)
  end if

  end subroutine
