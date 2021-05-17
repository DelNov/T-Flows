!==============================================================================!
  real function U_Tan(Flow, s)
!------------------------------------------------------------------------------!
!   Computes tangential velocity component at the cell near a wall.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  integer                   :: s     ! wall face
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c1, c2
  real                     :: u_tot_sq, u_nor_sq, u_nor
!==============================================================================!

  ! Take aliases
  grid => Flow % pnt_grid

  ! Take cells around the face
  ! (Note that c2 should be less than zero!
  !  Also note that there is no ckecking yet!)
  c1 = grid % faces_c(1, s)
  c2 = grid % faces_c(2, s)

  ! Compute tangential velocity component
  u_tot_sq = Flow % u % n(c1) * Flow % u % n(c1)  &
           + Flow % v % n(c1) * Flow % v % n(c1)  &
           + Flow % w % n(c1) * Flow % w % n(c1)
  u_nor  = ( Flow % u % n(c1) * grid % sx(s)     &
           + Flow % v % n(c1) * grid % sy(s)     &
           + Flow % w % n(c1) * grid % sz(s) )   &
         / grid % s(s)
  u_nor_sq = u_nor**2

  if( u_tot_sq  > u_nor_sq) then
    U_Tan = sqrt(u_tot_sq - u_nor_sq)
  else
    U_Tan = TINY
  end if

  end function
