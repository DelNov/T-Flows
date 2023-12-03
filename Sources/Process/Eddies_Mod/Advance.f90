!==============================================================================!
  subroutine Advance_Eddies(Eddies)
!------------------------------------------------------------------------------!
!>  The subroutine is designed to advance all eddies through an inlet plane in
!>  scale-resolving simulations.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * The subroutine first creates an alias, Flow, which points to the         !
!     pnt_flow field in the eddies structure. This field contains information  !
!     about the flow field in which the eddies are being advanced.             !
!   * It then iterates over all the eddies (n_eddies) in the eddies structure. !
!   * For each eddy, the subroutine updates its position (x, y, z) based on    !
!     its velocity (u, v, w) and the time step (dt) from the flow field. This  !
!     essentially moves the eddies through the plane in which it is defined.   !
!   * After updating the positions of the eddies, the subroutine checks if     !
!     any eddies have left the specified plane (x-plane, y-plane, or z-plane). !
!     This is done by comparing the distance of each eddy from the plane to    !
!     its length.                                                              !
!   * If an eddy has moved beyond its specified plane, it is replaced by       !
!     calling the Eddies % Place_Eddy subroutine. This ensures a continuous  !
!     generation of eddies in the domain, maintaining the turbulence           !
!     characteristics at the inlet.                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Eddies_Type), target :: Eddies  !! parent class; collection of eddies
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  integer                   :: e
!==============================================================================!

  ! Create an alias
  Flow => Eddies % pnt_flow

  ! Advance all eddies
  do e = 1, Eddies % n_eddies
    Eddies % eddy(e) % x = Eddies % eddy(e) % x  &
                         + Eddies % eddy(e) % u * Flow % dt
    Eddies % eddy(e) % y = Eddies % eddy(e) % y  &
                         + Eddies % eddy(e) % v * Flow % dt
    Eddies % eddy(e) % z = Eddies % eddy(e) % z  &
                         + Eddies % eddy(e) % w * Flow % dt
  end do

  ! Check which eddies left
  if(Eddies % x_plane < HUGE) then
    do e = 1, Eddies % n_eddies
      if( abs(Eddies % eddy(e) % x - Eddies % x_plane)  &
            > Eddies % eddy(e) % length) then
        call Eddies % Place_Eddy(e)
      end if
    end do
  else if(Eddies % y_plane < HUGE) then
    do e = 1, Eddies % n_eddies
      if( abs(Eddies % eddy(e) % y - Eddies % y_plane)  &
            > Eddies % eddy(e) % length) then
        call Eddies % Place_Eddy(e)
      end if
    end do
  else if(Eddies % z_plane < HUGE) then
    do e = 1, Eddies % n_eddies
      if( abs(Eddies % eddy(e) % z - Eddies % z_plane)  &
            > Eddies % eddy(e) % length) then
        call Eddies % Place_Eddy(e)
      end if
    end do
  end if

  end subroutine
