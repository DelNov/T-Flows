!==============================================================================!
  pure real function U_Tan(Flow, Grid, s)
!------------------------------------------------------------------------------!
!   Computes tangential velocity component at the cell near a wall.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), intent(in) :: Flow
  type(Grid_Type),   intent(in) :: Grid
  integer,           intent(in) :: s     ! wall face
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2
  real    :: u_tot_sq, u_nor_sq, u_nor
  real    :: u_wall_tot_sq, u_wall_nor_sq, u_wall_nor
!==============================================================================!

  ! Take cells around the face
  ! (Note that c2 should be less than zero!
  !  Also note that there is no ckecking yet!)
  c1 = Grid % faces_c(1, s)
  c2 = Grid % faces_c(2, s)

  ! Compute tangential velocity component
  u_tot_sq = Flow % u % n(c1) * Flow % u % n(c1)  &
           + Flow % v % n(c1) * Flow % v % n(c1)  &
           + Flow % w % n(c1) * Flow % w % n(c1)
  u_nor  = ( Flow % u % n(c1) * Grid % sx(s)      &
           + Flow % v % n(c1) * Grid % sy(s)      &
           + Flow % w % n(c1) * Grid % sz(s) )    &
         / Grid % s(s)
  u_nor_sq = u_nor**2

  u_wall_tot_sq = Flow % u % n(c2) * Flow % u % n(c2)  &
                + Flow % v % n(c2) * Flow % v % n(c2)  &
                + Flow % w % n(c2) * Flow % w % n(c2)
  u_wall_nor  = ( Flow % u % n(c2) * Grid % sx(s)      &
                + Flow % v % n(c2) * Grid % sy(s)      &
                + Flow % w % n(c2) * Grid % sz(s) )    &
                / Grid % s(s)
  u_wall_nor_sq = u_wall_nor**2

  if( u_tot_sq  > u_nor_sq) then
    U_Tan = abs(  sqrt(u_tot_sq - u_nor_sq)             &
                - sqrt(u_wall_tot_sq - u_wall_nor_sq))
  else
    U_Tan = TINY
  end if

  end function
